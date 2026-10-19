# Amended CUDA + CuPy development roadmap

Prepared 2026-10-19. This is an amendment, not a claim that the production
campaign or GPU fitter has been completed.

## 1. Location, authority, and purpose

This document is deliberately saved in the requested `LiouvilleCOPS-cuda`
directory. The existing implementation and original plan are in the different,
underscore-named directory [LiouvilleCOPS_cuda](../LiouvilleCOPS_cuda/).
Do not move, duplicate, or rename that implementation as part of this amendment.
Use its [original ROADMAP.md](../LiouvilleCOPS_cuda/ROADMAP.md) together with this
document; where they conflict, the changes below define the revised target.

The objective is to accelerate the restored fitting scheme on the SHAP5 788k
surface, including optional NG-only visibility and cliff activity enhancements,
without silently changing its physics, data selection, objective, or derivatives.
CUDA is an execution backend. It neither removes oscillations by itself nor
establishes that a recovered activity map is physically unique.

The principal amendment is to make the **free 29-parameter baseline and its
30/31-parameter physical-template extensions** the primary production target.
The fixed-baseline higher-degree harmonic experiments remain secondary modes.
The second essential amendment is to reject HEAL-corrupted gravity products and
admit only independently audited, correctly oriented gravity and cliff products.

Source records used for this amendment:

- [CUDA STATUS.md](../LiouvilleCOPS_cuda/STATUS.md),
  [VALIDATION.md](../LiouvilleCOPS_cuda/VALIDATION.md), and
  [GRAVITY_REFERENCE.md](../LiouvilleCOPS_cuda/GRAVITY_REFERENCE.md).
- [HEAL-corruption.md](../LiouvilleCOPS/HEAL-corruption.md) and
  [GRAVITY_REPAIR.md](../LiouvilleCOPS/GRAVITY_REPAIR.md).
- [RESTORATION_STATUS.md](../LiouvilleCOPS/RESTORATION_STATUS.md),
  [RESTORED_BASELINE.md](../LiouvilleCOPS/RESTORED_BASELINE.md), and
  [JOINT_TEMPLATE_FIT.md](../LiouvilleCOPS/JOINT_TEMPLATE_FIT.md).
- The current CPU evaluator, profile runner, gravity audit/cache, template
  selector, and their tests under `../LiouvilleCOPS/tools` and `tests`.

These records are a dated snapshot. Before implementation, inspect the live
status and artifact hashes; do not overwrite concurrent development or treat
historical “running” statements as current process status.

## 2. Why the original campaign emphasis must change

### 2.1 Restore parameter freedom, not just the old degree number

The useful archived result is in
`RES-COPS-101826-6pm-7459fd3c350733bb38c87799df62a3cdcbdf9b34`, including
`output_fit_full_current16`. Its fitted parameters were 25 degree-0–4 H2O
activity coefficients, three independent analytical temperature parameters,
and one CO2 multiplier. Lambda was zero. The selected data included saturated
NG and re-enabled late ExcludeNG observations.

The recent degree-4 controls instead had zero free parameters. The degree-5
experiments freed 11 higher-degree coefficients while fixing the original
baseline coefficients, temperature parameters, and CO2 multiplier. Those are
different inverse problems. A degree-4 fixed control is not a reproduction of
the archived 29-parameter fit, and the original baseline is not the archived
optimized baseline.

The archived fit's RMS log residuals were about 0.3121 NG and 0.1766 RG.
However, it reached the 1000-evaluation limit rather than converging, and its
temperature parameters included Tmax = 30 K and T75 = 400 K at their bounds,
with Tmin about 49.26 K. Thus the pressure agreement is a useful numerical
reference, not proof of a physically credible temperature solution. Preserve
these facts in comparison reports and do not present the archive as a validated
physical truth.

### 2.2 Do not identify lambda as the demonstrated cause of oscillations

Spatial activity variations project into time through changing illumination,
rotation, geometry, and instrument response. Strong spatial corrections can
therefore create temporal humps or dips. Changing which parameters are free can
matter much more than a small penalty on those parameters.

The recent fixed-baseline degree-5 scan did not provide a matched degree-5,
lambda-zero control establishing that regularization caused the oscillations.
Its curves changed little across the sampled positive lambdas. A selected
L-curve corner is not sufficient evidence that the selected fit is physically
preferred. The GPU port must not promise smooth curves, silently filter model
pressures, or change the optimizer to hide this behavior.

Use the restored lambda-zero baseline first. Investigate positive lambda,
higher harmonics, temperature ordering, and alternative selection profiles as
separately labeled controlled experiments.

### 2.3 Replace the invalid gravity reference

The HEAL defect used a nominal absolute perturbation too small to be reliably
represented at the metre-coordinate scale of this body. Self-intersection tests
could consequently misclassify facets and reverse already correct winding.
An offset-halving test can pass a consistently corrupted field; it is not an
orientation test.

The earlier 100k-versus-788k gravity rejection involved HEAL-affected gravity.
Retract its interpretation as a validated resolution comparison. Real geometry
differences remain, but their slope effect needs a corrected comparison. Keep
the exact 788k surface for the main campaign; a coarser gravity surface is not an
unannounced performance optimization.

## 3. Existing work to retain, and evidence still missing

| Area | Recorded state on 2026-10-19 | Amendment action |
|---|---|---|
| C01 environment | Isolated Python 3.11.14, CuPy 13.6.0, two RTX A6000 devices tested | Retain locked environment and diagnostics; do not reinstall merely because the original plan describes missing dependencies |
| C02 storage/backend | FP64 data layout, memory admission, identities, CPU facade validated | Extend configuration and parameter-state tests to the new joint 29–31 parameter modes |
| C03 gravity | Exact FP64 kernels, audited orientation, synthetic and 788k sample gates recorded as passed | Retain numerical evidence; bind it to exact code, mesh, and binary identities |
| C03 full field | `full_field.npz` recorded under `cuda_c03/gravity_20261019T171225Z` | Locate and verify artifact; it is not yet a production template bundle |
| CPU gravity repair | Topology audit, DISABLE path, cache invalidation, provenance guard implemented | Validate production cache import/rebuild and regenerated cliff templates |
| CPU joint fitting | Free baseline/T/kappa plus nonnegative templates implemented and tested | Make this the new GPU evaluator reference |
| CPU campaign preparation | Read-only B/N/C/NC command preparation implemented | Use it to freeze inputs, selections, and launch contracts |
| C00 acceptance | Tolerance policy/sign-off and production fit vectors remain pending in CUDA status | Reconcile with repaired CPU changes and obtain required acceptance; do not mark complete by inference |
| C06/C07 pressure fitting | CUDA status reports no pressure kernels yet | Implement and validate the revised primary modes before production claims |

CUDA status records a full-field median slope of about 23.1 degrees, area
fractions 13.56% at 45 degrees and 3.77% at 60 degrees, and all finite gravity
values. These are reported C03 results, not a new verification performed by this
document. Sample agreement does not constitute a full CPU/GPU field comparison.

The CPU restoration suite records 97 passed and seven skipped tests, with an
external eight-target C00 gravity check separately passed. Its synthetic C++
facet-reader replay is not a full mission SPICE replay. Preserve that distinction.

## 4. Revised scientific and numerical contract

### 4.1 Primary fit families

| Label | Free parameters | Purpose |
|---|---|---|
| B29 | 25 baseline coefficients + 3 temperatures + kappa | Restored baseline on the 788k mesh |
| N30 | B29 + one NG-only amplitude | Test the late visibility-region enhancement |
| C30 | B29 + one cliff amplitude | Test the slope-associated enhancement |
| NC31 | B29 + both amplitudes | Test their joint contribution and degeneracy |

Use degree 4 and enhancement degree 4 for these families. Do not pass the
legacy-activity, temperature, or kappa freezing options. Template amplitudes
are nonnegative; baseline coefficients retain their original freedom.
Keep the original degree-0–4 basis convention, ordering, and normalization.
Do not substitute the separately defined normalized higher-harmonic basis.

The current joint parameter order is baseline coefficients, selected templates,
temperatures, then kappa. Read names and scales from the CPU contract instead
of independently hardcoding an assumed column order in CUDA.

For H2O activity, evaluate

```text
u(f) = Y(f) a + A_scale B(f) c
A(f) = max(0, u(f))
```

There is ONE clip after addition. `A_scale` is the area-weighted RMS of the
original run's unclipped baseline map and remains fixed. Use the same derivative
branch as CPU: zero at and below the clip. Temperature and kappa derivatives
must correspond to the current joint state, not the original activity map.

The primary reproduction uses lambda zero, independent 30–400 K temperature
bounds, original initialization, and the existing TRF/scaling/tolerance settings.
The CPU preparation utility is the authority for executable commands. Extra
starts, ordered temperatures, or altered bounds are separate experiments.

### 4.2 Temperature, time dependence, and GPU reuse

Preserve the existing analytical temperature law, illumination/shadow branches,
floors, heliocentric scaling, species physics, and RG velocity correction.
Evaluate those terms and all requested derivatives at each new parameter vector.
Fixed-temperature precomputation remains valid only in explicitly fixed modes.
A cache must not reuse temperature-dependent factors merely because geometry
is unchanged. Include every free parameter and relevant mode in its identity.

The GPU port does not add thermal inertia, delayed sublimation, gas dynamics,
or a new gravity-driven production law. Nucleus rotation continues to enter
through the reference geometry/illumination and temperature/source calculations.
The cliff amplitude is an empirical test of enhanced activity on a prescribed
region, not a calibrated erosion or surface-renewal mechanism.

### 4.3 Data selection and objective

Preserve the exact profile flags, source rows, timestamps, pressure values,
saturation handling, exclusion rules, pressure floors, and gauge normalization.
Each gauge's log-pressure residuals are divided by the square root of that
gauge's selected sample count. Do not silently introduce censoring, weights,
interpolation, or late-time upweighting.

The recorded archive/current comparison has 930 NG and 1146 RG observations
in common. The 788k campaign includes four additional points per gauge,
yielding 934 NG and 1150 RG. Freeze explicit row identities: matching counts
alone does not prove matching data. Report common-observation comparisons
separately from the all-available production fit so a selection change is not
misattributed to mesh refinement.

### 4.4 Positive-lambda behavior, when explicitly requested

For joint modes, penalize the TOTAL unclipped change from the original prior:

```text
delta_u = Y (a - a_original) + A_scale B c
penalty = lambda * integral(delta_u^2 dS) / integral(u_original^2 dS)
```

Preserve baseline/template and template/template cross terms. Retain the CPU
PSD Gram-matrix square root, roundoff handling, and errors for materially
negative eigenvalues. Do not replace this with separate coefficient penalties
or add a ridge to conceal rank deficiency. Keep this small operation in CPU
FP64 initially. Starting from an archived optimum must not redefine the prior.

Warm-start support, resumable optimizer checkpoints, and eager per-lambda result
writing are not part of the first CPU restoration implementation. Implement and
test them explicitly if needed; do not describe them as existing functionality.

### 4.5 Gravity and mask identities

Audit the mesh before DISABLE evaluation: finite vertices, valid indices,
nonzero triangles, no duplicate faces, closed opposite-oriented edge incidence,
connected surface, manifold vertex links, Euler characteristic and positive
signed volume. The present audit does not prove absence of geometric
self-intersections; preserve trusted mesh provenance and document this limit.

Retain density 533 kg/m3, G = 6.67430e-11, self-gravity at centroid plus 0.1 m
outward normal, centrifugal acceleration at the centroid, and the reference
SPICE angular velocity/frame/epoch. Record the offset-halving audit separately.

Use the existing repaired cache algorithm identity:
`homogeneous_polyhedron_self_gravity_at_outward_offset_v2_disable_after_topology_audit`.
Use the actual orientation-policy constant from `gravity_topology.py` and
`integrity_check = DISABLE after topology audit`. A GPU importer must truthfully
record its own backend/compiler/evidence, not fabricate a CPU execution record.
If the schema cannot express that distinction, extend and test the schema.

The geometry-only NG mask is the union of `Kn > 0 AND KPhi == 0 AND KV == 0`
at the SAME point over 2016-09-30T10:33:10 through 10:39:26 UTC. Compute it from
original FP64 kernels. It is not “never seen by RG,” and does not represent a
source that is switched on only during that interval. Apply its fixed spatial
activity enhancement throughout the modeled trajectory.

Use the 60-degree cliff mask for the primary test and 45 degrees in a separate
sensitivity test. Record area, overlap, response coverage, and threshold-near
facets. Do not fit both nested cliff thresholds by default. NG masks require
geometry/interval/kernel identities; cliff masks require gravity/omega/threshold
identities. Neither is reusable solely by matching a mesh filename.

## 5. Amendments to the original stages

### C00 amendment — Re-freeze the repaired reference

1. Preserve the original fixtures as historical records. Add a versioned fixture
   set bound to the repaired CPU tree, dirty diff or commit, and exact inputs.
2. Add B29, N30, C30, and NC31 evaluator states: original start, perturbed interior
   states, zero/positive template amplitudes, clipping transitions, temperature
   branches, active bounds, and subsequently accepted production states.
3. Include selected-row tables, residual/Jacobian column names/scales, original
   prior, objective components, output exchanges, and configuration hashes.
4. Reconcile stale CUDA status statements about the running scan, HEAL builder,
   and pending CPU repair. Record status updates only after checking evidence.
5. Review the existing machine-readable tolerance policy. Explicitly approve
   pressure, log residual, Jacobian, objective, gravity-vector and slope budgets
   before tuning. Do not invent new numerical tolerances to obtain a pass.

Gate: reproducible repaired fixtures and an accepted policy. No corrupted cliff
fixture is promoted to truth; missing production states remain explicitly pending.

### C01 amendment — Retain the working environment

Retain and verify the existing environment lock and both device diagnostics.
Record runtime/driver/compiler/library versions in every production manifest.
Do not infer capabilities from the original pre-installation inventory.
The existing finding that cuSPARSE reductions are not bitwise repeatable
supports retaining custom fixed-order pressure/Jacobian reductions.

Gate: fresh diagnostics and supported CPU-reference interoperability, without
unnecessary dependency changes.

### C02 amendment — Extend contracts, not data physics

Add joint-mode state serialization, requested-column validation, bounds/scales,
template identities, and original-prior identity. Ensure all 31 columns can be
requested and returned in CPU order. Expand memory accounting to include joint
intermediates, derivative tiles, persistent state, and worker contexts.

Retain FP64 contribution storage, checked integer indices, int64 CSR offsets,
row ordering, and input corruption checks. GPU-memory estimates must describe
the live 31-column implementation, not only the existing static load experiment.

Gate: CPU facade equivalence for all four primary modes and rejection tests for
stale/mismatched templates, column sets, and parameter states.

### C03/C04 amendment — Admit corrected gravity into reusable products

1. Locate the recorded C03 full field and report through the configured output
   root. Verify hashes, shape ordering, dimensions, units, finite values,
   constants, audit, target positions, omega, and code/compiler identities.
2. Compare with corrected DISABLE CPU reference samples and difficult/threshold
   cases. Preserve the existing 832-target and synthetic evidence where valid.
3. Choose and document the production acceptance route: a corrected full CPU
   comparison, or an explicitly accepted exact GPU-field import supported by
   numerical/sample evidence. Sample validation is not full-field equivalence.
   If acceptance authority is unresolved, ask before declaring this gate closed.
4. Implement a checked importer/cache writer if reusing C03 output. A filename
   change or copied metadata is not a valid conversion.
5. Retain or implement two-GPU target partitioning and checkpoint/restart with
   atomic completion records, per-chunk checksums, identities, and no reuse of
   partially written chunks. Test interruption and mismatched restart refusal.
6. Keep HEAL-defect artifacts quarantined for provenance; never overwrite or
   quietly rehabilitate their cliff fields. Do not recompute exact gravity when
   a fully admitted cache already exists.

Gate: production-readable, correctly identified full-field cache plus documented
validation route. Two-GPU acceleration is not required merely to reuse a valid
one-GPU field.

### C05 amendment — Rebuild/admit templates independently of fitting

Regenerate effective slopes and cliff masks from the admitted gravity cache.
Retain the valid geometry-only NG mask, verifying its independent dependencies.
GPU Boolean reductions must be race-safe and exactly match CPU masks.

Verify threshold-sensitive facets against corrected FP64 CPU calculations and
record how ties are handled. The recorded C03 field had no facet within 1e-6
degrees of either threshold; verify this for the artifact actually imported.
Generate mask plots, area/overlap statistics, response-weighted diagnostics, and
a manifest linking every mask to its inputs. Response weights do not redefine
the Boolean source-region predicate.

Gate: the current CPU selector accepts repaired cliff provenance, all mask IDs
agree, and changing any dependent input invalidates reuse. An old bundle's
unaffected NG column may be selected without authorizing its cliff columns.

### C06 amendment — Prioritize the full free-parameter evaluator

Implement B29 first, then N30/C30/NC31 using the same fused or tiled evaluator.
Keep static geometric factors resident; recompute parameter-dependent physics.
Return pressures and analytic derivatives to existing CPU residual/regularization
and SciPy TRF code. Preserve residual-only versus Jacobian-request behavior.

Test all columns against CPU analytical derivatives and branch-aware finite
differences. Include zero-template equivalence, positive amplitudes, dependent
templates, gauge velocity terms, floors, clipped facets, and bound states.
Verify that free temperatures do not accidentally use fixed-temperature caches.

Gate: accepted FP64 budgets for all primary modes and unchanged C++ replay
tolerance. Small pressure error alone does not establish correct derivatives.
Higher-degree fixed-baseline modes may follow without blocking the primary modes.

### C07 amendment — Parallelize evaluations, not fictitious optimizer iterations

Keep one SciPy optimizer trajectory coordinated on the host. Use persistent
spawned GPU workers, split complete observations by contribution count/work, and
assemble small pressure/Jacobian outputs in original row order. Do not fork a
live CUDA context or assume the two devices form one automatic memory pool.

Benchmark one fit on two GPUs against independent fits on separate GPUs. Cap
memory and nested CPU pools explicitly. Use the requested 36 CPU threads for
appropriate parallel CPU work, including C++ replay; a serial optimizer does
not become a 36-process fit by changing OMP_NUM_THREADS. Gravity's TBB parallelism
must be controlled through its actual backend, not assumed to obey OpenMP.

Gate: one-/two-device objective equivalence, deterministic pinned reductions,
failure/cancellation cleanup, bounded memory, and measured latency/throughput.

### C08 amendment — Keep mixed precision optional

Retain FP64 as the release baseline. The archived inverse problem is ill
conditioned, clipped, and bound-active; A6000 FP32 throughput is not evidence
that FP32 is scientifically sufficient. Preserve FP64 branch decisions,
gravity, final reductions, residuals, derivatives as required by validation,
and small linear algebra. Test storage precision separately from arithmetic.

Gate: approved end-to-end budgets, including optimizer behavior and physical
conclusions, not only pointwise pressures. Failing this optional stage must not
delay a correct FP64 release or justify relaxing stopping tolerances silently.

### C09 amendment — Replace the primary production matrix

1. Run B29 with original initialization and lambda zero on the 788k data.
2. Run N30, C30, and NC31 with identical base settings and data; C/NC wait for
   repaired cliff admission. Preserve the original prior across all cases.
3. Compare common archived observations separately from all-available results.
4. Repeat important comparisons with the transition-cut profile. Exercise other
   supported profiles in integration tests without silently changing the main fit.
5. Run the 45-degree cliff sensitivity separately. Introduce positive lambda,
   ordered temperatures, additional starts, or higher harmonics only as labeled
   extensions with matched controls.
6. Produce full ten-hour and final-descent comparison plots, log residual plots,
   and source-region/activity diagnostics. Show excluded and saturated points
   distinctly rather than hiding them.
7. Report objective by gauge, RMS log errors, selected counts, stopping reason,
   nfev/njev, optimality, bound states, clipping, singular values/rank, parameter
   correlations, template overlap, and sensitivity to initialization/selection.
8. Replay selected final states through the full C++/SPICE forward model using
   fitted ActivityFit and facet enhancement exchange files. Verify the latter's
   hash binding to the fitted baseline, not the original prior.

Gate: numerical equivalence and a scientific comparison report. An optimizer
success flag, visually smooth curve, or better training residual alone is not
evidence that cliffs caused the observed NG jump. Quantify degeneracy and keep
the local-topography explanation a hypothesis unless supported independently.

### C10 amendment — Benchmark the revised problem honestly

Measure CPU, one-GPU and two-GPU implementations on identical B29 and NC31
states, data, precision, objective, and stopping criteria. Separate cold load/JIT,
warm pressure evaluation, all-column Jacobian evaluation, optimizer overhead,
gravity, masks, replay, and end-to-end wall time. Synchronize timings correctly.

Do not advertise an 11-column fixed-temperature speedup as the 31-parameter
free-temperature speedup. Distinguish latency from campaign throughput and
report host/device memory, contention, failures, and numerical deviations.
Retain CPU fallback and release only validated configurations.

## 6. Execution order and genuinely independent work

The following tracks can proceed concurrently once ownership and snapshots are
fixed; this is a development plan, not authorization to launch new production
jobs from this document.

- Reference/contracts: C00 updates and C02 joint-state fixtures.
- Gravity admission: C03 evidence reconciliation, C04 cache import/restart,
  then C05 cliff products.
- Evaluator: C06 B29 and NG/synthetic-template development using repaired CPU
  fixtures; it does not need a completed full cliff cache.
- Reporting/tests: B/N/C/NC manifests, diagnostics, output exchange tests,
  and performance harness without claiming production results.

Dependency sequence:

```text
repaired CPU reference + accepted numerical policy
    -> C02 joint contract -> C06 FP64 B/N/C/NC -> C07 multi-GPU
verified C03 field -> C04 production cache -> C05 repaired cliff bundle
    both numerical fitting and required template gates -> C09 production
    C09 + matched benchmarks -> C10 release
C08 mixed precision: optional branch after validated FP64, not a prerequisite
```

B29 and N30 production do not depend on gravity/cliff reconstruction. C30 and
NC31 do. Production GPU fitting requires evaluator gates even if CPU fits can
already run. Different workers must not concurrently edit a shared contract or
overwrite the same output/cache directory. Merge and rerun integration tests
at each dependency boundary.

## 7. Paths, artifacts, and reproducibility

Keep generated outputs outside both application source directories. Respect
the existing output resolver and `COPS_OUTPUT_ROOT`; the portable default is
the containing `AMPS/cops_output`, not a baked-in home-directory path.
The current machine's external output root may be selected explicitly; record
its resolved value without making it a universal installation requirement.

Current external inputs are under
`/data/vtenishe/rosetta_data_end_of_mission_cops`: `data_COPS_L2`, `sources`,
the shape models including `cg_spc_shap5_788k_cart.bdf`, and the persistent
`liouville_cops_polyhedral_gravity` package. These are deployment locations,
not proof of input identity. Retain command/config overrides and content hashes.

Each production output must include the full command, effective configuration,
code identity, original prior and start identities, free/fixed parameter list,
selection reasons, data/kernel/shape hashes, gravity/mask provenance, backend
precision and device identity, optimizer settings, stopping status, validation
results, and plot/statistics locations. Reusing caches must not reuse stale
optimizer outputs. Refuse accidental output overwrites.

## 8. Decisions and completion criteria

Do not fill these unresolved items with assumptions:

- Obtain any still-required C00 numerical-policy acceptance.
- Agree whether production gravity admission requires full corrected CPU-field
  comparison or accepts the validated exact GPU-field import route.
- If free temperature solutions remain at implausible bounds, request a separate
  physical-model decision; do not quietly impose ordering or a new thermal law.
- If new starts or checkpointing are required, specify their semantics and
  preserve the distinction between initialization and regularization prior.

The amended effort is complete only when the primary fit families are supported
with truthful provenance, repaired cliff products are reusable, FP64 CPU/GPU
evaluator and full C++ replay checks pass, the full-duration comparisons and
diagnostics exist, and measured performance is reported. This roadmap itself
does not launch fits, certify physical uniqueness, or mark pending stages done.
