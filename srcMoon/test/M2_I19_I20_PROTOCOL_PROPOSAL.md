# I19/I20 numerical verification protocol — detailed review proposal

**Approved by the user on 2026-11-13.** This document retains the reviewed
v3 design. The executable contract is [M2_I19_I20_PROTOCOL.json](M2_I19_I20_PROTOCOL.json);
implementation and measured evidence are tracked in STATUS.md. Approval does
not turn missing data/predecessors, short diagnostics or a planned cluster
allocation into formal qualification. Earlier proposal snapshots are retained
in the campaign artifacts. Scientific laws and old acceptance limits are not
relaxed by implementation.

Source basis: `06eafde4802f0344694612c8f3c21f2b28e1a09f`, plus the separately
recorded uncommitted D02 metadata/test/documentation work. Earlier versions
are retained in the review evidence directories. Version 3 incorporates the user-required
one-hour aggregate development limit and deferred cluster qualification. The
approved numerical design is consumed through the active protocol JSON; this
review document alone is not execution evidence.

## 1. Recommended decision package

| Decision | Recommendation | Reason / limit |
|---|---|---|
| Physical scope | Neutral Na and Ar, current literature coefficients, chemistry OFF; eight MPI ranks, OpenMP=1, individual weights ON. | Exercise the established M2 scope. Future species/conversions require their own coverage. |
| Clock | Common particle dt for both species, bounded by the smaller fully initialized production dt. Align the grid to physical checkpoints before running. | Equal iteration counts otherwise represent different species durations. |
| Refinement | Full 3x3 dt/weight design, with reference weights frozen independently of dt. | Separates integration sensitivity from statistical resolution and their interaction. |
| Physical cases | Source-fed zero-inventory case plus a separate controlled surface-cohort case. | The first tests normal startup/budgets; the second identifies residence and re-adsorption without censoring away the initial cohort. |
| Horizon | Short development horizon: 24 common coarse particle steps. Candidate cluster H = 2,550,240 s. | Development diagnostics cannot establish month-long retention; the long reporting interval is not an exact lunation or D02 averaging interval. |
| Replicates | Two explicitly selected development seeds for wiring/diagnostics; sixteen proposed cluster seed pairs per grid cell. | No sixteen-seed inference from two development replicates. The candidate cluster design has 288 runs. |
| New numerical margins | 1% of a frozen atom scale; one percentage point for fractions; 1% of H for restricted residence; 0.1% numerical-domain loss. | Engineering error allocations proposed for review, not published observational uncertainties. |
| Statistics | Paired-seed differences, simultaneous approximate 95% intervals over a predeclared family; confidence half-width at most half the margin. | A small or zero mean difference alone is insufficient. Sixteen seeds do not guarantee adequate precision. |
| Sites | At most 48 metadata-selected polar sites, before model evaluation. | Stratified verification coverage; not a representative global observational sample. |
| Execution | Whole development campaign <=3600 s, including builds/plots/cleanup. Long qualification is explicit cluster work of order a few hundred CPU-core hours. | User-required execution policy; full scientific closure remains gated. A proposed 300 core-hour cluster target needs a measured feasibility estimate. |

## 2. Current implementation and changes actually needed

The current supporting path `Moon::IntegrationTests::RunM2BudgetCheck` calls
ordinary `amps_time_step`, sources, force mover, surface callbacks and
`Moon::M2Budget`. The ledger independently sums owned gas particles and a
replicated surface reservoir counted once, including pending surface fluxes.
`Exosphere::SourceProcesses::InjectFiniteSurfaceInventory` uses integrated
competing thermal/PSD hazards and actual partial-packet weights.

Relevant current code:

```cpp
// srcMoon/main_lib.cpp, after the initial species initialization messages
PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
    _AR_SPEC_,_NA_SPEC_,0.1);

// srcMoon/MoonIntegrationTests.cpp, current 20-step supporting fixture
const double dt=std::min(saved[0],saved[1]);
PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]=dt;
PIC::ParticleWeightTimeStep::GlobalTimeStep[_AR_SPEC_]=dt;
rnd_seed(731993+PIC::ThisThread);
```

Thus a campaign seed in a Python manifest is insufficient: this hard-coded
seed must be replaced by an explicit, verified runtime selection in the new
formal path. Existing supporting/reference modes must retain their declared
behavior. Likewise, use final initialized timesteps, not the earlier printed
initialization messages. The retained I19 receipt reports
`0.6855865301040126 s`, despite earlier messages showing approximately
`6.855865 s` before the Ar copy.

`MoonM2Budget.cpp` currently writes completed steps, species, counts and atom
terms; it needs physical epoch/elapsed time and window endpoints. The evolving
external environment in `main_lib.cpp` advances with species zero's dt.
The dedicated common-clock mode must verify gas motion, surface injection,
thermal advancement and SPICE/frame advancement use the same physical interval.
Changing only one array or merely relabeling a step column is insufficient.
Preserve the normal production dispatch; do not replace it with a second test
integrator.

## 3. Common time, weights and initial conditions

Let dt_bound be the minimum Na/Ar production timestep after initialization,
including all species copies. For eight equally spaced reporting intervals,
set

```
m = ceil(H / (8 * dt_bound))
N0 = 8 * m
dt0 = H / N0
```

Run refinements dt0, dt0/2 and dt0/4. The step counts are N0, 2*N0 and 4*N0;
checkpoints H*j/8 are integer step indices in every run. Compute timestamps
from start epoch plus integer-indexed elapsed time, with the documented
floating-point error bound, rather than repeatedly summing a large TDB epoch.
Record the actual dt and duration. Do not interpolate inventories or assume
an unsupported shortened final mover step. Validate local block/source dt
views after installing the global control.

A planning example using the retained 0.6855865301040126 s bound gives
N0 = 3,719,800, dt0 about 0.685585246519 s and **14,879,200 finest steps**.
This arithmetic is not a forecast that the present fixture can run that long.
The actual formal dt_bound must be measured for its frozen input.
At that illustrative bound the whole two-case design contains
672*N0 = 2,499,705,600 completed eight-rank timesteps. Weight-dependent particle
work is additional to this step count. Do not remove the legacy Ar 0.1 factor
solely to meet a runtime target; any alternative timestep restriction needs
its own physical/numerical justification and review.

For each species, preserve W_start from ordinary baseline initialization at
full precision. A separate resource/precision pilot may select a finer
W_ref = W_start / 2^b, with nonnegative integer b. Freeze b and W_ref before
formal execution; choose them for numerical precision/capacity, never model-
observation agreement. The formal weight axis is W_ref, W_ref/2, W_ref/4.
Do not rerun the source-particles-per-step weight heuristic when dt changes.
Source intensities, SI inventories and physical time remain fixed.

Freeze all geometry, AMR/domain extents, source coefficients/switches, force
settings, epoch/frame conventions, thermal initial state and forcing across
the grid. Use a separately qualified thermal solution with an unchanged
thermal accuracy setting; retain the approved I17 error budget. Particle
refinement cannot replace thermal, terrain or AMR convergence.

Two declared physical cases are recommended:

1. **Source-fed startup:** preserve the approved zero initial surface stock,
   the existing normal gas-initialization policy and all input-declared Na
   sources; Ar has the continuous uniform 3.8e20 atoms/s source. Inventory is
   transient. Do not claim equilibrium or M1 closure from this case.
2. **Controlled surface cohorts:** zero initial gas, external injection OFF,
   explicit Na/Ar stock at the frozen sites at the common t0. Recommend
   N_site,s = 32*W_start,s, extending the existing controlled fixture scale.
   Freeze these physical atom counts before selecting W_ref, and preserve
   them for every dt/weight run. Existing enabled finite-stock release laws
   stay active; switch-specific thermal-only controls belong in separate
   analytical fixtures. This is a verification initial condition, not an
   observation or claimed lunar abundance. Released atoms can return, stick
   and release again through production paths.

Both cases use the same actual initial thermal state/forcing and horizon.
Require explicit evolving-environment mode for these coupled physical cases;
retain frozen forcing only in separately labelled analytical controls. The
thermal initial state must be qualified at the actual t0. Do not reuse a
June-ended spin-up state at an unrelated January epoch or force an epoch
choice from D02's unresolved averaging metadata.
No monthly forcing history or averaging epoch is inferred from D02. Production
thermal initialization must pass its declared gates before formal coupling.

## 4. Conservation gates: keep the existing requirements

For every completed ledger window, and separately for Na and Ar, require:

- Exact **gas model-particle count** closure. A continuous surface stock has
  no model-particle count that can be added to the gas count.
- Gas, effective surface, combined atom and kilogram closure within the
  existing declared production-arithmetic bounds. Independently reconstruct
  the identities from serialized terms; do not fit bounds to residuals.
- Converted-in/out exactly zero in the chemistry-OFF case.
- Separate external injection, signed prescribed reservoir supply, adsorption,
  credited deposition, thermal/PSD release, physical escape and numerical loss.
- The existing strictly-below-0.5% dominant-source imbalance gate for its
  nonzero-source windows. A closed cohort with no external source still has
  exact conservation; do not divide by zero or invent a source.

Internal adsorption/release cancels in the combined atom identity. Uncredited
adsorption must have an explicit supported physical-sink interpretation and
its own entry; it cannot be silently called numerical loss or ignored.
The I09 energy classifier remains the declared lunar two-body exit classifier,
not a claim of ultimate three-body escape.

The existing supporting analyzer's denominator takes a maximum that also
includes internal transfer totals. Preserve its reviewed check, expose that
scope, and add an explicitly named external-source normalization to formal
source-fed reporting. Large gross recycling must not conceal an imbalance
relative to external supply. Freeze the detailed formal equation/denominator
in I19's contract before implementation; a proposal is not permission to
reinterpret the accepted 0.5% criterion silently.

Check the arithmetic bound's usefulness over long windows. If cumulative
operation counts exhaust or make it uninformative, fix/verify the summation
or window accounting; do not widen a physical acceptance threshold. Preserve
separate cumulative and completed-window columns.

## 5. Proposed additional numerical margins

Define a common, per-species physical normalization

```
A_s = initial gas atoms + initial surface atoms
      + integral of the predeclared expected external supply over [0,H].
```

Freeze this value and source integral before formal refinement; do not use
a different random realized supply or model residual as each run's divisor.
Validate the source integral independently to an error negligible relative to
the proposed budget. For the source-free cohort case A_s is its initial stock.
An exactly empty/source-free control uses exact zero invariants.

| Quantity | Proposed comparison criterion |
|---|---|
| Gas/surface atoms and integrated injection, adsorption, release and escape | Each normalized refinement difference within +/-0.01 of A_s. |
| First-release survival and current surface occupancy | Each absolute fraction difference within +/-0.01. |
| Restricted mean first residence | Difference in mu_H/H within +/-0.01; about 25,502.4 s for this H. |
| Numerical-domain loss | Each run individually loses at most 0.001*A_s; keep numerical/physical losses separate. This is a proposed new long-case ceiling, not a relaxation of the current supporting fixture's zero-loss rule. |
| Statistical precision | Simultaneous confidence half-width <=0.5 of the relevant margin, as well as the complete interval lying inside the margin. |

One percent is a proposed resolution budget. It is not an observational
uncertainty, a proof of an infinite-resolution error bound or evidence that a
small source process is accurate to 1% of its own flux. Report small transfer
terms in absolute units as well, and verify each source/rate against its own
analytical law in controlled cases. Do not let a large external-supply divisor
make an unexercised path look verified. The original 0.5% conservation gate and
roundoff closure remain stricter bookkeeping obligations.
Also retain gas/surface changes relative to their own populations. A passing
A_s-scaled gate alone is not a claim of 1% relative gas-density accuracy when
most atoms are stored elsewhere. Science-specific relative-error requirements
must be frozen for the corresponding observables before observational use.

## 6. Refinement and statistical design

This formal statistical design applies to the deferred cluster tier. The
short two-seed development cases have their own diagnostic scope in section 9.

For each of the two cases, use all nine dt/weight combinations. Use sixteen
seed indices j=0..15 with base seed `731993 + 1009*j`; for the MPI-only profile,
propose rank seed = base + rank. Verify complete reseeding of the configured
uniform and cached derived-distribution state. Retain the RNG mode, effective
rank seeds and initialization receipts; `src/general/rnd.{h,cpp}` is the
production authority. Distinct seed integers alone do not prove independent
streams. Audit the active generator and seed/stream policy before describing
replicates as statistically independent; never replace the core RNG silently.

Pair seed indices across refinements. No bitwise identity, guaranteed common
random-number coupling or inverse-square-root weight scaling is assumed.
Particles, timesteps, facets and correlated checkpoints are not independent
replicates.

Freeze these primary comparisons before evaluating them:

- All 12 adjacent axis pairs in the 3x3 grid, plus each of the eight non-finest
  cells compared directly with the finest cell; remove the two duplicate
  pairs. This gives **18 unique grid-pair comparisons**.
- For both cases, six global metrics (gas, surface, integrated external
  injection, adsorption, release, escape), two species and eight nonzero
  checkpoints: at most 3456 scalar comparisons across the 18 pairs.
- For the cohort case, three site endpoints at H (first-release survival,
  current surface occupancy, mu_H/H), two species and at most 48 sites:
  at most 5184 additional comparisons.
- Reserve g_max = **8640** for the formal numerical confidence family. Freeze
  the actual site/metric registry before AMPS results; keep this conservative
  divisor even if an observational stratum is unavailable. Missing required
  strata remain explicit coverage gaps, not numerical PASS.

For normalized paired differences d_j, recommend

```
mean(d) +/- t_(1 - 0.05/(2*g_max), 15) * sd(d)/sqrt(16).
```

Require the whole interval inside the equivalence margin and its half-width
within the precision budget. This applies paired-difference inference and
Bonferroni simultaneous coverage to a fixed comparison family
([NIST paired observations](https://www.itl.nist.gov/div898/handbook/prc/section3/prc311.htm),
[NIST Bonferroni method](https://www.itl.nist.gov/div898/handbook/prc/section4/prc473.htm)).
The simultaneous confidence interpretation is conditional on valid marginal
intervals; it does not establish normality or stream independence.

At g_max=8640 and 16 seeds the critical t is approximately 6.82061. To obtain
a normalized half-width of 0.005, the paired-seed standard deviation must be
at most about 0.00293229. This planning arithmetic illustrates why sixteen
seeds may be insufficient; the 20-step supporting fixture does not establish
the required full-duration precision.

**Statistical qualification is required.** Sixteen replicates with Student-t
are an approximation for non-normal simulation outputs; finite variance alone
does not guarantee 95% coverage. Sparse, highly skewed or unsupported all-zero
stochastic results cannot receive zero-error PASS from this formula. A
normality-test non-rejection is not proof. Require a documented justification
from aggregate-event support and an independent statistical preflight; if
unsupported, preserve FAIL / insufficient statistical qualification and propose
a separately frozen larger-seed or justified alternative design. Do not
silently switch methods, append selected seeds or widen margins after a miss
([NIST discussion of mean confidence limits](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/conflimi.htm)).

Analytically verified deterministic identities/structural zeros use their
exact oracle, not a stochastic zero-variance shortcut. Keep every replicate;
no cherry-picking. Store all 18 comparison intervals, confidence assumptions,
coverage diagnostics and explicit fail reasons.

For the timestep axis, reject a resolved increase in sensitivity on refinement;
allow an unresolved change only when both comparison intervals already meet
the fixed numerical and precision budgets. For the weight axis, require
resolution equivalence and precision; do not demand monotone means from
individual random realizations. No global second-order convergence assumption
is made merely because one mover is second order.

## 7. Residence, re-adsorption and warming: distinct physical quantities

Use passive mass tags on the controlled initial cohorts. Record actual release
and surface-transfer events without an extra RNG draw. Preserve the physical
mixed-reservoir representation: tags carry fractions of real packet/stock
atoms, and total tagged atoms never exceed actual atoms. Re-adsorption does not
reset an atom's first-release history. Verify tag conservation independently.
Do not create a parallel surrogate mover or inject fictitious uncounted atoms.

For a cohort initially at a site at t0, report:

- S_first(t): fraction that has **never left** its initial adsorption episode.
- F_surface(t): fraction presently on any surface, including returned atoms;
  optionally retain original-site occupancy as a separate diagnostic.
- First-release fraction 1-S_first(t), distinct from the gross accumulated
  release flux, which can exceed initial stock after repeated cycles.
- mu_H = integral from 0 to H of S_first(t) dt, the restricted mean first
  residence; atoms still resident at H remain in this integral.
- Escape, numerical loss and current gas fractions, with a closing tagged
  atom identity. Later residence episodes are separate statistics.

The endpoint cohort quantities have the proposed simultaneous refinement
limits. Save their intermediate curves as diagnostics without treating the
checkpoints as extra independent replicates or claiming an unfrozen endpoint
as a primary PASS.

At constant temperature and thermal release alone, the independent reference
is S_first(t)=exp(-t/tau), and mu_H=tau*(1-exp(-H/tau)), evaluated stably with
expm1. The no-release limit is mu_H=H. This follows the exponential survival
law ([NIST exponential distribution](https://www.itl.nist.gov/div898/handbook/eda/section3/eda3667.htm)).
For changing temperature or additional finite-stock release channels, integrate
lambda_total(t), with S_first(t)=exp(-integral(lambda_total)dt), using an
independent converged quadrature; do not use a mean temperature inside the
Arrhenius exponential. Freeze the independent quadrature tolerance before
comparison and demonstrate its error <1% of the applicable numerical margin.
This verifies kinetics conditional on the supplied temperature history; it
does not independently qualify the thermal model that produced that history.
Retain current U12 exact/roundoff criteria for its constant-hazard checks.

Extend production-path controlled fixtures to cover constant cold/warm states,
piecewise warming with an exactly known hazard integral, and a smooth warming
history with independently converged quadrature. Choose temperatures by the
unchanged species laws and declared dimensionless H/tau targets, not by
observational agreement. These are imposed analytical verification forcings,
not lunar observations. Keep thermal-only and combined-channel checks explicit.
Censoring is information, not a reason to omit residents or substitute a mean
computed only from releases
([Stan survival/censoring documentation](https://mc-stan.org/docs/stan-users-guide/survival.html)).

The existing tau>=H flag is only a thermal-timescale diagnostic. With constant
thermal hazard, tau=H implies survival exp(-1), about **36.8%**, not permanent
or even 90% retention. PSD/other active loss hazards can reduce actual survival
further. Keep geometric PSR, thermal-timescale flag and actual kinetic
retention separate. Do not introduce an irreversible sink or claim infinite
residence from a finite-horizon result.

## 8. Polar sites and external-data limits

Recommend 2 hemispheres x 4 absolute-latitude bands (65-75, 75-80, 80-85,
85-90 degrees) x PSR/non-PSR x 3 sites, at most 48. Use half-open bands, with
90 degrees included in the last band. Require native-valid D02/D03 geometry,
explicit frame/datum conversions and the existing eight-neighbour ambiguity
rule. Select the lowest SHA-256 values of an unambiguous canonical UTF-8
product-ID/line/sample tuple within each stratum. Recommend a compact JSON
array `[product_id,line_1_based,sample_1_based]`, UTF-8, no whitespace, no trailing
newline; break a digest tie by product ID and then integer line/sample. Publish
that exact encoding,
selected native records and hashes before model evaluation. Co-registration
must meet a predeclared coordinate/footprint rule, not an undocumented nearest
value selected after comparison.

Missing strata or native gaps stay unavailable. If multiple selected products
map to one production facet, retain the mapping and flag the shared model
support; do not call them independent model sites. The three sites/stratum
are coverage probes, not proof that the whole stratum is resolved.

D03 PRP annual-mean temperature at 2 cm and annual-maximum surface temperature
remain separate observables. Neither is an epoch-resolved surface temperature
history. No invented heating curve, uncertainty or kinetic-temperature
conversion is supplied from those summaries. D02 remains PARTIAL until its
exact temporal/geometry contract and I15 reference are resolved. Metadata-only
site preparation and supporting analytical fixtures may proceed; formal I20
must wait for READY exact inputs and I15/I17/I18/I19 plus same-campaign I13.
The existing June/December selection and hold-out protections remain unchanged.

## 9. User-directed execution tiers and resource policy

The user requires executable development tests within one hour and permits
longer qualification to be deferred to a larger cluster for a few hundred
CPU hours. These are execution requirements. The quantitative physics/
statistics proposals elsewhere remain unapproved. No new executable flag or
aggregate watchdog is claimed to exist already.

### Development: one aggregate hour

Budget **3600 s end to end**, not per test. Target 3000 s of planned work,
leaving 600 s for variability and controlled shutdown/reporting. Include required
configuration/clean compilation, local/linked tests, data QA, plots, manifests,
report aggregation and cleanup. A proposed planning allocation is:

| Work | Target seconds |
|---|---:|
| Input/dependency/data preflight | 180 |
| Required clean configuration/build variants | 1500 |
| Local kernels, baseline/production-path checks, short refinement diagnostics | 960 |
| CSV/PNG/EPS and manifest aggregation | 180 |
| Planned cleanup/report finalization | 180 |
| Reserved headroom/termination | 600 |
| **Whole command deadline** | **3600** |

These allocations are targets to measure on the development host, not proof
that the current code meets them. Prefer an isolated checkout and known cached
external dependencies, so a deadline cannot leave the user's application
reconfigured. Time source/config/header hashes and compilation separately;
reuse a binary only with verified matching source/build identity. Do not hide
necessary compilation outside the hour to claim an end-to-end timing PASS.

Development content should retain U/local analytical checks, required M0
regressions, read-only D01/D02/D03 QA, normal eight-rank terrain/surface/Na-Ar
budget paths, and seed/clock/weight/cohort observer checks as implemented.
Preserve the short supporting fixture; extend it with a proposed common short
horizon H_dev=24*dt0. This allows all eight checkpoints to land on integer
steps at dt0, dt0/2, dt0/4: **24/48/96 steps**. All weight refinements compare
the same physical interval. Two disjoint development seed indices from
2000003+1009*j, j=0,1, give a candidate two-case 3x3 design of **36 short runs**,
with up to two concurrent eight-rank jobs after measuring aggregate memory.
Measure this workload before freezing it as the default timed suite.

The short sweep evaluates production wiring, exact signed closure and declared
analytical fixtures. Its resolution/seed differences are diagnostics, not the
formal sixteen-seed confidence gate or evidence of month-long equilibrium/
retention. The 500-step science sampling length remains unchanged; a 24-96-step
case has no completed 500-step science sample. Existing completed-window
storage fixtures test that separate path. Controlled warming analytical cases
retain their explicit verification forcings and do not stand in for lunar
observations or long thermal spin-up.

Use the existing global framework and I19/I20 supporting namespaces, without
new permanent IDs. Every short case retains its scope, runtime, seed, input/
source/build/data hashes, closing ledger and comparison files. Do not omit a
selected case after seeing results. Qualify the complete runtime on the host;
if selected work cannot finish within the budget, preserve ERROR/incomplete
outputs and revise the execution design explicitly before a new attempt.

The aggregate deadline must include builds and all child processes. Before
launching a phase, pass its **remaining** budget rather than resetting a
3600-s timer. Start controlled termination at the reserved deadline, preserve
completed outputs, reap the full MPI/build process groups and write the final
summary within the hard limit. Cleanup must not start an unbounded restore
build. The current runner only passes `test_runtime_seconds` separately to
I14/I18/I19; this does **not** implement the required campaign-wide guard and
must be fixed before advertising a guaranteed one-hour command.

The earlier v2 plan for twenty planning jobs of 300 s each excluded builds and
grace. It is superseded. Any extra capacity/precision pilot is part of the same
hour, or a separately selected cluster task; it is never silently appended to
the development command.

### Cluster qualification: explicit selection and measured core-hours

Defer long physical-duration transport, thermal spin-up, full dt/weight/seed
convergence and formal retention to this tier. Keep the original science
horizons, parameters and criteria unless a distinct scientific revision is
reviewed. A missing data/operator/predecessor remains SKIPPED. Cluster access
or development PASS does not satisfy it.

Recommend **300 CPU-core hours as an aggregate planning target**, adjustable
before launch within the user-authorized few-hundred-hour scale. This is not
300 hours for every seed, process or retry. Use

```
allocated_core_hours = sum(allocated_cores * wall_seconds) / 3600.
```

For eight MPI ranks and one thread/rank, one job-hour consumes eight allocated
core-hours. Thus 300 core-hours corresponds to 37.5 summed eight-rank job-hours,
or ideally 18.75 wall hours with two such jobs continuously active. Record
scheduler-reported actual CPU and charged allocation separately, since their
definitions can differ. Include failed attempts and retries; do not reset a
campaign's accounting by splitting it into job arrays.

The full proposed two-case 3x3/16-seed design remains **288 long simulations**,
not a default development workload. At the illustrative retained clock it has
about 2.50 billion timesteps. A 300 core-hour allowance would permit only about
7.81 minutes on average per eight-rank simulation, including its allocated
setup/postprocessing. No existing measurement establishes that feasibility.
A larger machine alone does not guarantee lower core-hour cost.

Pilot the physical cases and smallest-weight/highest-occupancy settings;
measure compilation, initialization, thermal work, post-initialization
throughput, particle growth, load balance, allocated cores and output size.
Require a forecast within the selected total allowance before scheduling.
If it exceeds a few hundred core-hours, present a separately reviewed staged
qualification design or keep the expensive remainder deferred. Do not shorten
physical time, discard seeds, loosen precision or grant PASS to an incomplete
matrix merely to fit the allocation. Formal-run batching must preserve the
same protocol/source/executable/data identities and all attempts.

### Common outputs, resource controls and statuses

Retain the proposed aggregate development RSS ceiling **200,000,000,000 bytes**
and artifact ceiling **50,000,000,000 bytes**, with named host memory/disk
reserves. Begin with one job until measured allocation permits two. Cluster
memory/disk limits must be explicit in its submission manifest, not copied
blindly from the workstation. Production allocation failures retain the
normal `exit(__LINE__,__FILE__,msg)` behavior with detailed remediation.

Keep particle/surface restart writes and full-volume mesh output OFF; retain
500-step sampling where applicable, physical checkpoint/closing ledgers,
resource histories and CSV/PNG/EPS comparisons after each completed stage.
Avoid exhaustive particle/event trajectory dumps. M1's 5000-step rule is
unchanged; short M2 diagnostics do not complete M1.

Report development and cluster qualification separately:

- PASS: the declared selected scope was evaluated and met its own criteria.
- FAIL: a completed selected case missed a physical/numerical/statistical gate.
- ERROR: a selected case/build suffered corruption, infrastructure/resource
  failure or campaign-deadline termination; preserve finished case results.
- SKIPPED: the cluster tier was not selected, or a required capability, data
  product or predecessor is unavailable. Deferred is never qualification PASS.

The last short I19 run used about 21.26 GiB sampled peak RSS and 23.73 s for
20 steps including setup, under
`test_output/srcMoon/m2-next-20261113T042634Z/m2-development/m2-development/supporting/I19/`.
It informs development planning but does not predict month-long precision,
host-wide end-to-end timing or cluster core-hours. A current development test
result must state measured whole-command runtime and whether compilation was
performed, reused with verified identity, or unavailable.

## 10. Implementation sequence and review boundary

After scientific review, implement in this order:

1. Update the existing I19/I20 contracts and source-roadmap amendment
   consistently; retain the external roadmap read-only and document both
   hashes/authority in ROADMAP_RECONCILIATION.md.
2. Add explicit protocol/runtime controls following `src/INPUT_FILE.md`:
   case/seed selection, common dt, pinned per-species atom weights, horizon/
   checkpoints, cohort stock/site table, output/resource settings. Snapshot
   all resolved values in `.in` and protocol receipts. Names/syntax need the
   normal parser implementation; this proposal does not claim they exist.
3. Extend the production observer/ledger with physical time, signed terms,
   window/cumulative checks and verified passive cohort accounting. Preserve
   baseline modes and no-extra-RNG diagnostic behavior.
4. Add independently checked constant/warming references, clock/seed/weight/
   rank-ownership protocol fixtures and the audited source/mover path matrix.
   Every advertised active path must be exercised in a controlled positive
   case. Physically rare events in the nominal case must retain an independent
   rate/oracle check; a zero count is not rate verification. Finite-inventory
   sputtering remains explicitly unsupported until its distinct law is
   specified; chemistry conversions and hybrid OpenMP are outside this scope.
5. Implement the fixed-design replicate launcher and independent scorer;
   register them through existing I19/I20/global entries. Preserve every
   attempt, fail reason, seed, input/source/build/data hash and resource trace.
6. Implement and verify the aggregate development deadline, then run bounded
   supporting checks within the whole one-hour envelope. Freeze actual reference
   weights, registry, statistical assumptions and measured cluster cost before
   formal launch.
   Rerun M0 and the required regressions after any shared production change.
7. Run formal qualification only when all numerical/data/operator predecessors
   are satisfied. No I13 reference is installed automatically, no December
   residual is inspected and no observational agreement is inferred from
   conservation or analytical verification.

The user-directed execution amendment changes the source roadmap SHA-256 to
`f1a2f9265f97f36da9ef9f6a17b3c83783d9f4b51882e6d81081ef7d567757cd`;
external roadmap SHA-256 remains
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source amendments govern the approved supporting-development exception;
common scientific gates govern qualification; native labels and external
AGENTS govern data. Only the user-directed execution-budget amendment is
adopted in the source roadmap; scientific criteria are unchanged.

**Approval boundary:** the user approved the design and requested implementation.
No further approval for these same choices is required. Long execution still
requires measured cost, genuine formal predecessors, exact READY inputs and
frozen statistical/site/normalization evidence; no missing factual metadata
is supplied by that approval. See STATUS.md for actual scope and results.
