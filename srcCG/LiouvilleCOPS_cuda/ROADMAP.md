# LiouvilleCOPS CUDA + CuPy development roadmap, revision 2

| Revision | Date | Content |
|---|---|---|
| 1 | 2026-10-19 | Original plan, kept unchanged as `ROADMAP_v1.md` |
| 2 | 2026-10-19 | Incorporates `AMMENDED_ROADMAP.md` (the revised fitting scheme and gravity admission) and the review comments in `AMENDMENT_REVIEW.md` |

Where this revision and `ROADMAP_v1.md` differ, this revision governs. The amendment
remains the source document for its decisions. Stage evidence is recorded in `STATUS.md`;
the numerical policy is in `VALIDATION.md`.

**Stage status at this revision:**

| Stage | Status |
|---|---|
| C00 | complete; C00b re-freeze complete |
| C01 | complete |
| C02 | complete, including the joint-mode extensions |
| C03 | complete (full-field admission passed in C04) |
| C04 | single-GPU complete: 100k and 788k equivalence passed, admitted caches written; two-GPU in C07 |
| C05 | complete: visibility masks and restart bundle (bitwise equal to the CPU builder) |
| C06 | complete: numerical gate passed on the C00b subset and the full 1150-point set; Jacobian budget approved |
| C07 | step 1 done: unchanged fitter on one GPU (first B29 fit, 17.5 min); two-GPU work next |
| C08 – C10 | not started |

## 1. Objective and scope

Port the expensive gravity, masking, pressure and analytic-Jacobian work of the
self-contained `LiouvilleCOPS` application to CUDA, orchestrated from Python through CuPy.
Keep SciPy's bounded trust-region least-squares optimizer, and keep the CPU application as
the definition of every scientific quantity.

**Primary target (amended):** the *free* 29-parameter baseline fit and its 30/31-parameter
physical-template extensions on the SHAP5 788k surface (§3.3).
* The fixed-baseline higher-degree harmonic experiments (degrees 5–8) are secondary,
  labelled modes. They are not the campaign's centre.
* Only independently audited, correctly oriented gravity and cliff products are admitted.
  HEAL-derived products are rejected (§3.1).

CUDA is an execution backend. It does not remove oscillations by itself, and it does not
make a recovered activity map physically unique. Any change in fitted-curve behaviour must
come from the fitting scheme, documented and tested against the CPU reference, never from
GPU-side smoothing.

Dataset (counts are data, never program constants):
* 787,754 facets;
* 1,150 observation points and 245,515,533 point–facet contributions;
* all-available selection: 934 NG and 1,150 RG observations, of which 930 NG and 1,146 RG
  are common with the archived fit.

Initial transport-kernel construction, ray tracing and SPICE stay on the CPU.

## 2. Environment (C01, as built)

* **Python and packages:** Python 3.11.14 in `.venv` with CuPy 13.6.0; NumPy 1.24.4 and
  SciPy 1.10.1 (identical to the CPU reference); polyhedral-gravity 3.3. Pinned with hashes
  in `requirements-cuda.lock`.
* **CUDA:** system toolkit 12.4 with driver 550.
* **GPUs:** two RTX A6000 (cc 8.6, NVLink NV4). Both report the kernel-execution watchdog
  as enabled, so long kernels are split into bounded launches.
* **Reference interpreter:** the CPU reference keeps running under `/usr/bin/python3`
  3.8.10.
* **Diagnostics:** `tools/cuda_doctor.py`.
* **Reductions:** cuSPARSE reductions are not bitwise repeatable, so pressure and Jacobian
  reductions use custom fixed-order kernels.

## 3. Scientific behaviour that must be preserved

### 3.1 Gravity and slope

* **Evaluated quantity:** homogeneous self-gravity of the exact 788k polyhedron in metres;
  density 533 kg/m³; G = 6.67430e-11 m³ kg⁻¹ s⁻².
  * Self-gravity is evaluated at `centroid + 0.1 m · outward normal`.
  * The centrifugal term −ω×(ω×r) is evaluated at the centroid.
  * Slope = acos(−n̂·g_eff/|g_eff|).
* **Angular velocity:** ω is the SPICE body-frame angular velocity of the reference, with
  its epochs and identities.
* **Orientation:**
  * **No HEAL.** The package's ray test flips about 46% of the faces of a kilometre-scale
    mesh in metres (`GRAVITY_REFERENCE.md` §5.1).
  * A body is admitted only by the independent audit:
    * the CPU's `gravity_topology.require_outward_closed_body`: finite, valid indices,
      unique non-degenerate triangles, closed opposite-oriented edges, one connected
      component, manifold vertex links, Euler characteristic 2, positive volume;
    * the CUDA-side `gravity_body` audit as a cross-check.
  * It is then evaluated with `DISABLE`. The audit does not prove absence of geometric
    self-intersection, so mesh provenance remains required.
* **Cache identity:** algorithm
  `homogeneous_polyhedron_self_gravity_at_outward_offset_v2_disable_after_topology_audit`,
  with orientation policy `DISABLE after independent outward closed manifold audit v1`. A
  GPU-produced field records its own backend, compiler and evidence; it never copies CPU
  execution metadata (§5 C04).
* **Cliff templates:** binary thresholds ≥ 60° (primary) and ≥ 45° (separate
  sensitivity). Never fit both nested thresholds together by default. Density and
  thresholds are prescribed hypotheses, not fitted parameters.
* **Physical audits, not numerical allowances:** the offset-halving audit (512 samples,
  0.25° gate) and any resolution comparison. An offset-halving pass does *not* prove
  orientation: a consistently corrupted field passes it.
* **The earlier 100k-vs-788k "rejection"** used HEAL-corrupted gravity. It is retracted as
  a resolution study; the exact 788k surface remains the gravity boundary.
* **HEAL products** stay quarantined (`*_heal_defect`) and are never rehabilitated.

### 3.2 Visibility and masks

* `late_ng_not_rg` is the union over the selected points of the closed interval
  2016-09-30T10:33:10 – 10:39:26 UTC of facets with `Kn > 0` and, *at the same point*,
  `KPhi == 0` and `KV == 0`.
  * Computed from the original FP64 kernels with race-safe integer OR.
  * It is not "never seen by RG", and it is not a source switched on during the interval.
    Its fixed spatial enhancement applies over the whole modelled trajectory.
* **`late_rg_not_ng`** (added by the user's decision of 2026-10-19) is the RG counterpart:
  facets that only RG can see and NG cannot. It is the union over the same selected points
  of facets with `KPhi != 0 or KV != 0` and, *at the same point*, `Kn == 0`. These are the
  CPU builder's own per-point predicates (`point_rg = KPhi != 0 | KV != 0`,
  `point_ng = Kn > 0`), used in the opposite direction.
  * Interval: the same declared late interval by default, configurable.
  * Use: it may enter a fit as a nonnegative template, like `late_ng_not_rg`.
  * Construction: same FP64 kernels, same race-safe OR, same identities.
* Response-weighted maps are diagnostics. They never redefine the Boolean predicate.
* Reuse requires identities:
  * NG masks: geometry, interval and kernel identities;
  * cliff masks: gravity, ω, threshold and definition identities.

  A matching mesh filename is never enough.

### 3.3 Fitting scheme (amended)

`../LiouvilleCOPS/tools/fit_cops.py` remains the definition of the evaluator: species
sensitivities, the production law, shadow and SZA handling, heliocentric scaling, the
analytical temperature branches, gas masses, NG conversion, each RG pressure model, the
spacecraft-velocity correction and CO2 scaling.

**Primary fit families** (degree 4, enhancement degree 4, no `--fix-*` options):

| Label | Free parameters | Templates | Depends on gravity |
|---|---|---|---|
| B29 | 25 degree-0–4 H2O coefficients + Tmax, T75, Tmin + κ | none | no |
| N30 | B29 + 1 amplitude | `late_ng_not_rg` | no |
| C30 | B29 + 1 amplitude | `cliff_ge_60deg` | yes: repaired C05 bundle |
| NC31 | B29 + 2 amplitudes | both | yes: repaired C05 bundle |

**Labelled extension (user decision of 2026-10-19):** `late_rg_not_ng` (§3.2) may be
added as a further nonnegative template, alone or with the others, using the same joint
machinery: one clip, joint regularization, no frozen parameters. It is reported as a
labelled family (for example R30 = B29 + RG-only, NR31 = B29 + NG-only + RG-only), never
silently merged into the primary matrix. Its overlap and amplitude correlation with the
other templates are reported.

**Activity field.**

```text
u(f) = Y(f)·a + A_scale · B(f)·c        A(f) = max(0, u(f))
```

* There is exactly one clip, applied after addition. The derivative is zero at and below
  the clip (same branch as the CPU).
* `A_scale` is the area-weighted RMS of the original run's *unclipped* baseline map, and it
  is fixed.
* Template amplitudes are dimensionless and ≥ 0. Baseline coefficients are unbounded.
* The degree-0–4 basis keeps its original convention and normalization. The separately
  normalized higher-harmonic basis is not substituted.
* Temperature and κ derivatives use the current joint field.

**Parameter order and scaling** are read from the CPU `Problem` (names, `x0`, bounds,
scales), never hard-coded. The current order is: baseline coefficients, templates,
temperatures, κ. The C02 contract returns columns in exactly that order.

**Primary settings:**
* λ = 0;
* independent temperature bounds 30–400 K;
* the original initialization (run coefficients and temperatures, zero amplitudes);
* the existing TRF settings: `x_scale='jac'`, ftol, xtol and gtol 1e-10, max_nfev 1000.

The CPU preparation tool `tools/prepare_restored_campaign.py` is the authority for the
executable commands. Extra starts, ordered temperatures and altered bounds are separate,
labelled experiments.

**Data and objective:**
* Preserved exactly: profile flags, source rows, timestamps, pressures, saturation
  handling, exclusions, pressure floors.
* Log-pressure residuals of each gauge are divided by √(selected count of that gauge).
* Not introduced: censoring, extra weights, interpolation or late-time upweighting.
* Row identities are frozen explicitly. Matching counts do not prove matching data.
* Common-observation comparisons with the archive are reported separately from the
  all-available production fit.

**Positive λ (only when explicitly requested):** the penalty is the area integral of the
squared *total* unclipped change from the original prior,
`Y(a − a_orig) + A_scale·B·c`, divided by `∫u_orig² dS`.
* It includes the baseline/template and template/template cross terms, through the CPU's
  eigen square root of the Gram matrix.
* No added ridge.
* It is computed on the CPU in FP64.
* The prior is never redefined by a warm start.

**Temperature dependence:** temperatures are free in all primary modes. Every
temperature-dependent factor and derivative is evaluated at each new parameter vector.
Fixed-temperature precomputation is valid only in modes where the temperatures are
explicitly fixed. Cache keys include every free parameter and mode.

**Not added by the port:** thermal inertia, delayed sublimation, gas dynamics, or a
gravity-driven production law. The cliff amplitude is an empirical test of enhanced
activity on a prescribed region.

**Secondary modes:** the fixed-baseline degree-5–8 harmonics, positive-λ scans, ordered
temperatures, alternative profiles and multiple starts. They are supported or refused
explicitly and never block the primary modes.

### 3.4 Oscillations of the fitted curves

**Cause.**
* Spatial activity structure projects into time through rotation, illumination, viewing
  geometry and instrument response.
* With the baseline, temperatures and κ frozen, every misfit, including smooth trends, can
  only be absorbed by short-wavelength harmonics. Those appear as temporal humps and dips.
* The primary families limit spatial freedom to degree ≤ 4 plus at most two fixed-shape,
  nonnegative templates.
* λ was not shown to be the cause: the completed L5 scan had no matched λ = 0 control.

**What the port must not do:**
* promise smooth curves;
* filter model pressures;
* change the optimizer to hide behaviour.

**Diagnostics** (reported for every C09 fit and for the archived B29; reported only, never
a gate, as decided by the user on 2026-10-19). Per gauge, on the selected observations in time order:
* sign changes and lag-1 autocorrelation of the log residuals;
* second-difference roughness energy of the model log pressure against that of the data,
  over common times, for the full 10 h and the final descent separately;
* the clipped-area fraction along the optimizer trajectory.

**Credibility trigger.**
* The archived B29 ended at max_nfev with Tmax = 30 K, T75 = 400 K and Tmin ≈ 49 K. That
  is an inverted temperature law, which can trade off against the activity map.
* If the 788k B29 again ends at a temperature bound, C09 reports it and runs the
  ordered-temperature mode as a labelled control.
* It then stops for the user's physical decision. No ordering or new thermal law is
  imposed silently.

## 4. Architecture (built and planned)

```text
LiouvilleCOPS_cuda/
  ROADMAP.md (this file), ROADMAP_v1.md, AMMENDED_ROADMAP.md, AMENDMENT_REVIEW.md
  README.md, RUN.md, STATUS.md, VALIDATION.md, GRAVITY_REFERENCE.md, validation/tolerances.json
  cops_cuda/
    built:   errors, provenance, reference, formats, precision, memory, backend,
             cuda_sources, device_store, manifest, environment, gravity_body, gravity (GPL)
    kernels: store_audit.cu, gravity.cu (GPL)
    planned: gravity_cache_io (C04), masks + kernels/masks.cu (C05),
             evaluator + kernels/pressure.cu, jacobian.cu (C06), device_pool / workers (C07)
  tools/
    built:   c00_capture_reference, cuda_doctor, c02_device_load, c03_gravity_check
    planned: c00b capture, c04 cache import/build, c05 masks, c06 evaluator check,
             benchmark_cuda, fit runner integration (C09)
  tests/
```

* **Shared CPU code:** readers, profiles, `Problem`, regularization and reporting are
  imported from the bound reference tree (`cops_cuda/reference.py`), never copied.
* **Device use:** CuPy owns devices, arrays, streams and NVRTC compilation. Fused indexed
  reductions are RawModule kernels.
* **Multi-GPU:** persistent spawned workers, one per GPU; never fork after CUDA has been
  initialized.
* **Outputs:** written below `$COPS_OUTPUT_ROOT` in unique directories. The default
  resolver root is `<AMPS>/cops_output`; this host overrides it to
  `/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output`.

## 5. Stages and acceptance gates

Each stage delivers:
* commented code, focused tests and documented commands;
* machine-readable evidence;
* a `STATUS.md` entry.

Sample agreement is never reported as full-field equivalence, and pending production
evidence is never inferred.

### C00: freeze the reference. COMPLETE; C00b PENDING

* **Done:** fixtures for identities, masks, fit and gravity (`c00_20261019T064536Z`). The
  numerical policy was approved on 2026-10-19.
* **C00b (amendment): re-freeze against the repaired CPU tree.** It changed 13:02–13:06 on
  2026-10-19, after the C00 capture.
  1. Keep the C00 fixtures as history. Add a versioned set bound to the repaired tree's
     digest and exact inputs.
  2. **States:** add B29/N30/C30/NC31 evaluator states:
     * the original start and perturbed interior states;
     * zero and positive amplitudes;
     * clipping transitions and temperature-floor and lit/night branches;
     * active bounds;
     * the accepted production states as they become available.

     C30/NC31 states use synthetic templates until the C05 bundle exists. No corrupted
     cliff column becomes a fixture.
  3. **Context:** record the selected-row tables, parameter names and scales, the original
     prior, objective components, exchange outputs and configuration hashes.
  4. **Policy:** extend the approved budgets to these states. New tolerances need the
     user's approval; none are invented to obtain a pass.
* **Gate:** a reproducible repaired fixture set.

### C01: environment. COMPLETE

Keep the lock and the doctor. Record runtime, driver, compiler and library versions in
every manifest. Change dependencies only through the lock.

### C02: contract, layout, provenance. COMPLETE; joint extensions PENDING

* **Done:**
  * the contract: `ParameterState`, `ColumnRequest`, `EvaluationResult`,
    `ReferenceBackend`, and a `ModelFacade` that runs the unmodified CPU `Problem`;
  * the SoA FP64 store with int64 offsets;
  * memory admission with 20% headroom;
  * identities.
* **The NC31 joint mode already passes:** 31 columns, joint eigen-regularization, bitwise
  through the facade (`tests/test_backend_contract.py`).
* **Remaining:**
  1. Bind template and original-prior identities in `backend.identity()`, and refuse a
     stale or mismatched template, column set or parameter state.
  2. Serialize parameter states and requests, for C00b fixtures and restart records.
  3. Memory accounting for the *live* 31-column evaluator: joint intermediates, derivative
     tiles, persistent state, worker contexts. This is delivered with C06.
* **Gate:** CPU-facade equivalence for B29/N30/C30/NC31, plus the rejection tests.

### C03: exact FP64 gravity on one GPU. Numerical gate COMPLETE

* **Done:**
  * the GPLv3 translation;
  * precompute bitwise equal to the package;
  * ≤ 0.33 eps·M on the synthetic probes; ≤ 2.3e-13 relative on 832 + 32 targets at 788k;
  * 160 difficult points; NaN parity;
  * bitwise determinism across chunking and both GPUs;
  * the full 788k field: 3150 s, 13.56% of area ≥ 45°, 3.77% ≥ 60°, no facet within 1e-6°
    of a threshold (`cuda_c03/gravity_20261019T171225Z`).
* **Remaining (amendment):** bind the evidence to exact code, mesh and binary identities
  (the report manifest has them). Also use the reference audit `require_outward_closed_body`
  in `gravity_body`, keeping the local audit as a cross-check.

### C04: production gravity cache (amended)

1. **Verify the C03 artifact:**
   * hashes, facet order, units and finite values;
   * constants, audit, target positions, ω, and code and compiler identities.
2. **Acceptance route (decided by the user on 2026-10-19).** The gravity model is tested
   at full-field scale with the *repaired CPU builder* on the smaller-resolution nucleus
   `SHAP5_Volume.bdf`. That is its 100,000-triangle Gmsh surface, wound inward, so it is
   used through a documented exact global winding reversal (§8). The 788k transport facets
   are the targets: 787,754 targets × 100k faces.
   * The GPU evaluates the same body and targets.
   * The comparison is facet by facet under the approved budgets: ≤ 1e-12 relative,
     ≤ 1e-9° slope, threshold band re-resolved on the CPU.
   * This proves the evaluator on every target. The production 788k field then enters
     through the GPU import route, supported by that proof plus the 788k sample evidence
     (832 + 160 targets).
   * A by-product: the corrected 100k-vs-788k slope comparison that replaces the retracted
     HEAL-era rejection.
3. **Cache writer:** use the CPU's own `gravity_cache.write_cache` to write a
   repaired-algorithm cache.
   * The definition (algorithm, geometry identity, shape hash, density, offset) is
     backend-independent.
   * The calculation metadata records the backend, compiler, device and evidence, plus the
     two fields the unchanged CPU builder checks on reuse: `orientation_policy` and
     `topology`, equal to `gravity_topology.require_outward_closed_body`'s result, which
     `GravityBody.reference_audit` keeps verbatim. It also records the offset-halving
     audit.
   * The builder checks no package binary identity, so **no change to `../LiouvilleCOPS`
     is needed**. Changes there must not alter its workflow (user constraint, 2026-10-19).
   * Never copy CPU execution metadata.
4. **Two-GPU work:**
   * target partitioning by measured chunk time, keeping the original output order;
   * atomic per-chunk completion records with checksums and identities;
   * resume only verified chunks of the same configuration;
   * tests for interrupted writes, corrupt chunks, mismatched restarts, moved inputs,
     changed ω (derived fields recomputed) and reordered shapes (refused).
5. **Reuse and restart:** never recompute exact gravity when a fully admitted cache
   exists. Two-GPU acceleration is not needed merely to reuse a valid field. A restarted
   campaign or fit reloads the admitted gravity cache and the C05 mask bundle after
   identity checks (user decision 2026-10-19). Optimizer-state checkpoints come later.
* **Gate:** a production-readable, correctly identified full-field cache with a documented
  validation route; the CPU loaders accept it.

### C05: masks and effective fields (amended)

1. **Cliff masks:** regenerate effective slopes and cliff masks from the admitted cache.
   * Ties: a facet whose slope equals the threshold counts as ≥ the threshold.
   * Facets within 1e-6° of a threshold are recomputed on the CPU.
   * The threshold-band count is recorded for the imported artifact.
2. **NG and RG masks:** keep the geometry-only NG mask after verifying its own
   dependencies. It may be selected from the quarantined bundle (as the CPU tool allows);
   its cliff columns may not. Add the **RG-only mask** `late_rg_not_ng` (§3.2), computed
   from the same kernels and interval. The CPU builder does not produce it, so C05 is its
   producer. A CPU NumPy replica of the predicate is the test oracle, and it must agree
   bit for bit.
3. **GPU predicates:** race-safe integer reductions, bit-identical to the CPU masks.
4. **Bundle (the restart product):** written with the CPU's own
   `activity_templates.write_bundle`, holding `late_ng_not_rg`, `late_rg_not_ng`,
   `cliff_ge_45deg` and `cliff_ge_60deg`. It comes with plots, area and overlap statistics,
   response-weighted diagnostics, and a manifest linking every mask to its inputs.
   * The bundle metadata records the repaired gravity definition: algorithm `…_v2…`, the
     orientation policy, and `integrity_check = "DISABLE after topology audit"`. The
     unchanged CPU fitter's selection gate therefore accepts the cliff columns, while
     `backend` and evidence fields state truthfully that the masks were made by the CUDA
     backend.
   * No change to `../LiouvilleCOPS` is needed (user constraint, 2026-10-19).
   * Overlap of `late_ng_not_rg` with `cliff_ge_60deg` is recorded, for NC31
     identifiability.
* **Gate:**
  * the CPU selector accepts the repaired cliff provenance;
  * all mask identities agree;
  * changing any dependency invalidates reuse.

### C06: single-GPU FP64 evaluator (amended priority)

1. **Order:** **B29 first**, then N30/C30/NC31 with the same evaluator. The fixed-baseline
   harmonic modes follow and do not block.
2. **Physics:** keep the static geometric factors resident and recompute all
   parameter-dependent physics per call, including u = T^-½ and du/dT, because
   temperatures are free. A fixed-temperature specialization is allowed only for
   explicitly fixed modes, with keys that include the mode.
3. **Jacobian:**
   * one fused pass per column tile (`memory.plan_tiles`);
   * per-point fixed-order segmented partial sums, with multi-block rows reduced in a fixed
     order;
   * the masked factor `[u>0]·w·u` is shared by the baseline and template columns;
   * never a contribution × column array.
4. **Interface:** return pressures and requested columns only. Residuals, regularization
   and TRF stay in the CPU code through the C02 facade. Test both the residual-only and the
   Jacobian paths.
5. **Tests:**
   * every column against the CPU analytic derivatives;
   * branch-aware central differences;
   * zero-template equivalence, positive and dependent templates;
   * velocity terms, floors, clipped facets and facets crossing the clip, bound states;
   * proof that free temperatures never hit a fixed-temperature cache.
* **Gate:**
  * the approved FP64 budgets (pressure ≤ 512 eps·M against the CPU, derived log-residual
    and objective bounds);
  * the C++ replay tolerance (1e-6) unchanged;
  * the Jacobian budget, decided here from per-entry magnitudes and approved by the user.

### C07: two-GPU fitting (amended)

1. **One optimizer:** one SciPy trajectory on the host. Persistent spawned workers own
   devices and static data. Observations are split by contribution count and measured
   work; small results are assembled in the original row order.
2. **Benchmarks:** one fit on two GPUs versus independent fits per GPU, with capped memory
   and nested CPU pools.
3. **CPU threads:** 36 threads only for genuinely parallel CPU work, such as C++ replay. A
   serial optimizer is not parallelized by thread settings. TBB is controlled through its
   own backend (the affinity mask).
* **Gate:**
  * one-device and two-device objective equivalence;
  * deterministic pinned reductions;
  * clean failure and cancellation;
  * bounded memory;
  * measured latency and throughput.

### C08: mixed precision (optional)

FP64 is the release baseline. The inverse problem is ill-conditioned, clipped and has
active bounds, so FP32 throughput is not evidence of scientific sufficiency. Storage
precision is tested separately from arithmetic precision. A budget needs the user's
approval, and failing this stage never delays a correct FP64 release.

### C09: production matrix (replaced)

1. **Runs:**
   * B29 with the original initialization and λ = 0;
   * then N30, C30 and NC31 with identical settings and data;
   * C/NC only after C05;
   * the original prior is kept in every case.
2. **Comparisons:**
   * common archived observations (930 NG / 1146 RG), reported separately from the
     all-available results (934 / 1150);
   * repeat key comparisons with the transition-cut profile;
   * a 45° cliff sensitivity run on its own.
3. **Labelled extensions:** positive λ, ordered temperatures, extra starts and higher
   harmonics, each with matched controls.
4. **Report:**
   * objective by gauge, RMS log errors, counts;
   * stopping reason, nfev and njev, optimality;
   * bound states and clipping;
   * singular values, rank and correlations;
   * template overlap and amplitude correlation;
   * initialization and selection sensitivity;
   * **the §3.4 oscillation diagnostics**;
   * full 10-hour and final-descent plots with excluded and saturated points shown.
5. **Temperature bounds:** apply the §3.4 trigger if any fit ends at a temperature bound.
6. **Replay:** replay selected final states through the full C++/SPICE forward model with
   the fitted `ActivityFit.txt` and the facet enhancement file. Verify its hash binding to
   the *fitted* baseline.
* **Gate:** numerical equivalence plus a scientific comparison report. A success flag, a
  smooth curve or a lower training residual alone is not evidence that cliffs caused the
  NG jump.

### C10: performance and release

* **Measure:** CPU, one-GPU and two-GPU on identical B29 and NC31 states, data, precision,
  objective and stopping criteria.
* **Separate:** cold load and JIT, warm pressure evaluation, all-column Jacobian, optimizer
  overhead, gravity, masks, replay and wall time, with correct synchronization.
* **Do not mix cases:** never advertise a fixed-temperature or 11-column speedup as the
  31-parameter free-temperature speedup.
* **Report:** memory, contention, failures and numerical deviations.
* **Release:** keep the CPU fallback and release only validated configurations.

## 6. Numerical acceptance policy

`VALIDATION.md` / `validation/tolerances.json`, approved 2026-10-19:
* exact IDs and masks;
* pressure ≤ 64 eps·M (exact sum) and ≤ 512 eps·M (CPU);
* derived log-residual and objective bounds;
* gravity ≤ 16 eps·M, plus ≤ 1e-12 relative on campaign targets, plus NaN parity;
* slope ≤ 1e-9°; threshold band 1e-6°;
* bitwise repeatability;
* no mixed-precision budget.

**Approved at the end of C06 (2026-10-19):** the Jacobian budget, ≤ 1024 eps·M_J per entry.
**Provisional:** the end-to-end fit criteria (decided after the first full CPU-vs-GPU fit). New modes and states (C00b) use the
same budgets. Changes need a documented, user-approved amendment.

## 7. Execution order

The tracks below can proceed concurrently once ownership is fixed. This plan authorizes no
production run by itself.

* **Contracts and evaluator:** C00b fixtures → C02 joint extensions → C06 B29 → C06
  N30/C30/NC31 → C07. B29 and N30 need no gravity.
* **Gravity admission:** C03 reconciliation → [decision: acceptance route] → repaired full
  CPU field and comparison (or accepted import) → C04 cache → C05 repaired cliff bundle.
* **CPU campaign** (the user's, outside the CUDA stages): B and N fits can run on the CPU at
  any time with the commands from `prepare_restored_campaign.py`. C and NC wait for C05.

```text
repaired CPU tree + approved policy ─► C00b ─► C02 ext ─► C06 (B29 → N/C/NC) ─► C07 ─┐
C03 field ─► [acceptance route] ─► C04 cache ─► C05 cliff bundle ────────────────────┼─► C09 ─► C10
C08 (optional, after validated FP64)                                                  ┘
```

Different workers never edit a shared contract or write the same output or cache
directory concurrently. Integration tests are rerun at each dependency boundary.

## 8. Decisions

**Decided by the user on 2026-10-19:**
* **Gravity acceptance:** full-field CPU-vs-GPU test on `SHAP5_Volume.bdf` (100k faces)
  with the repaired CPU builder; the 788k field is then imported (C04 item 2).
* **Derived surface:** `SHAP5_Volume.bdf`'s 100k-triangle surface is consistently inward.
  It is used through an exact global winding reversal, written by
  `tools/derive_outward_surface.py` to
  `/data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/derived/SHAP5_Volume_surface_outward.bdf`
  with a provenance record. Vertices are byte-identical; nodes 2 and 3 are swapped in every
  face.
* **Oscillation diagnostics (§3.4):** reported only, never a gate.
* **Restart:** restart reuses persisted products: the admitted gravity cache (C04) and the
  NG/RG/cliff mask bundle (C05), each identity-checked on load. Saving and resuming the
  optimizer state (to adjust and restart a fit) is deferred to a later stage.
* **CPU tree constraint:** additions to `../LiouvilleCOPS` must not change its workflow,
  and substantial code changes are made only in `LiouvilleCOPS_cuda`. The C04 cache and
  C05 bundle are written with the CPU's own writers and accepted by its unchanged loaders.
* **RG mask:** `late_rg_not_ng`, facets only RG can see and NG cannot (§3.2), over the
  same late interval as `late_ng_not_rg`. It is part of the restart products and may be
  used as a fitting template (§3.3, labelled).
* **Order of work:** C00b and C02 extensions, then C06 B29. C05 masks and the C04 writer
  run alongside.

**Still pending (not to be filled by assumption):**
1. The physical decision if B29 again ends at temperature bounds (§3.4 trigger).
2. The end-to-end fit budget (after the first CPU-vs-GPU fit). The Jacobian budget was
   approved on 2026-10-19 (≤ 1024 eps·M_J per entry).

The amended effort is complete when:
* the primary families are supported with truthful provenance;
* repaired cliff products are reusable;
* the FP64 CPU/GPU evaluator checks and the full C++ replay pass;
* the full-duration comparisons and diagnostics exist;
* performance has been measured and reported.

## 9. Reference sources

* **CPU evaluator:** `../LiouvilleCOPS/tools/fit_cops.py` (`Model.prepare`, `Model.evaluate`,
  `Problem`, `fit_one`).
* **Gravity:** `tools/gravity_topology.py`, `tools/gravity_cache.py`,
  `tools/build_physical_activity_templates.py`.
* **Templates and campaign preparation:** `tools/activity_templates.py`,
  `tools/prepare_restored_campaign.py`.
* **CPU records:** `../LiouvilleCOPS/RESTORATION_STATUS.md`, `RESTORED_BASELINE.md`,
  `JOINT_TEMPLATE_FIT.md`, `GRAVITY_REPAIR.md`, `HEAL-corruption.md`, `HARMONIC_BASIS.md`,
  `ACTIVITY_ENHANCEMENT_FILE.md`.
* **This tree:** `AMMENDED_ROADMAP.md`, `AMENDMENT_REVIEW.md`, `GRAVITY_REFERENCE.md`,
  `VALIDATION.md`, `STATUS.md`, `ROADMAP_v1.md`.
* **External:**
  * CuPy kernels: https://docs.cupy.dev/en/stable/user_guide/kernel.html
  * polyhedral-gravity: https://esa.github.io/polyhedral-gravity-model/ (GPLv3; see
    `GRAVITY_REFERENCE.md` §7)
