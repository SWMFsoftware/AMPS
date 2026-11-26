# Proposal to complete lunar chemistry qualification

**Recommendation.** Complete I21 first using the commissioned sodium producer.
In parallel, implement the atomic part of I22 and qualify a narrowly defined
I23 helium pilot. Spend one bounded source-review session on the preferred
charge-exchange channel. Treat the complete molecular family, the approved
charge-exchange scope and same-campaign M1/M2 closure as three independent
dependencies for formal M3 completion.

This proposal follows section 5.43 of [STATUS.md](../STATUS.md#543-i21-homogeneous-execution-and-i22i23-preparation--2026-11-24-utc).
It builds on the existing [M3 implementation proposal](M3_COMPLETION_PROPOSAL.md),
[approved statistical and numerical addendum](M3_QUALIFICATION_ADDENDUM_20261118.md)
and [targeted physical-input audit](M3_PHYSICAL_INPUT_FOLLOWUP_20261120.md).
The source basis is `267fdc25fba6f9b510143425f390958801a64cd1` with the
local work recorded in section 5.43. Proposed case selections require review
and a frozen plan before execution.

## 1. Starting evidence and formal predecessors

| Gate | Evidence available | Outstanding requirement |
| --- | --- | --- |
| I21 Na photoionization | Six completed native launches across three time steps and two MPI layouts. Eighteen supporting checks passed, one strict numerical check failed, and 6,048 P02 checks passed. The 43 reactions per realization share the same 128 original parent histories. | Resolve the numerical failure and execute the complete source, shadow, mover, weight, MPI and statistical matrix. |
| I22 multi-species photon chemistry | D04 rates and branches, shared conservative event machinery, and 70 passing atomic reference-state checks for ⁴He, ²⁰Ne and ⁴⁰Ar. | Admit source/state-compatible atomic events and all seven H₂O and four OH branches, including charged breakup. Complete the U15 and U24 gates. |
| I23 electron impact | D05 Voronov coefficients, rate and analytic log-slope verification, and a shared reduced reservoir evaluator. U16 has 302 passing checks. | Bind a bounded source/state-specific provider to real particle events, then qualify the He, Na, Ne and Ar cases in the frozen protocol. |
| I24 charge exchange | Verified rate readers and candidate partial and total-capture sources. Shared conservative two-body and forward-approximation services are available. | Admit a source-supported product mapping and angular law or error bound, then verify actual collisions and ledgers under U17 and I24. |

The formal I21–I24 results are **four SKIPPED**. I11 failed independence in
the retained bounded M1 qualification, I12 passed, and I13 awaits a passing
I11 fine-grid case and review of the reference it produces. Complete I13,
then qualify I14–I20 in the same campaign, including R2.6 conservation. For
I15/I20, seek a timestamped independent illumination benchmark with
producer-defined ephemeris, solar-disc, Earth-occultation and observer-geometry
conventions. If these cannot be obtained for the frozen D02 products, bring
any alternative evidence and selection rules to roadmap review before
inspecting comparison residuals. Preserve the frozen products and their
PARTIAL status. For I17, keep the approved independent December spin-up,
frozen hold-out and numerical convergence rule. Run the remaining I19/I20
qualification at its planned resources.
[M1 completion checklist](../STATUS.md#71-exact-m1-completion-checklist)
[M2 remaining work](../STATUS.md#72-m2-work-remaining-in-addition-to-m1)

The [source roadmap](../development_and_validation.md) SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The canonical data-root roadmap at
`/data/vtenishe/moon_validation_data/development_and_validation.md` is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source amendments govern work order and resource tiers, including the
approved early M3 development exception and the prescribed Maxwellian pilot.
The shared scientific criteria, data definitions, test IDs and predecessor
requirements govern qualification. [Roadmap reconciliation](ROADMAP_RECONCILIATION.md)

## 2. I21 — Finish the existing sodium production campaign

1. **Preserve and diagnose the constant-hazard result.** The maximum
   normalized errors at 8, 16 and 32 steps were approximately
   `7.38e-17, 7.38e-17, 9.09e-17` in both layouts. Native hexadecimal event
   times show the same plateau before CSV projection. Preserve that FAIL and
   use the cohort to diagnose exponential-clock and ledger behavior. Record
   the arithmetic floor, each tested invariant and where rounding first
   enters. Apply the existing strict rule to the retained case. Any change in
   its role in the formal case design requires a separate review of the
   contract and this evidence.
2. **Add a varying-hazard case with an independent clock reference.**
   Recommend an illumination-uninterrupted, zero-force radial path with
   `r(t) = r0 + u*t`, `r0 > 0`, `u > 0` and source rate
   `k(t) = k1AU*(AU/r(t))²`. Its cumulative hazard is
   `H(t) = k1AU*AU²*t/(r0*(r0 + u*t))`. Invert `H(t) = -log(U)` using
   independent high-precision arithmetic and check event position/time,
   survivor classification and the daughter’s unused step. Use a prescribed
   analytical Sun-distance fixture to exercise the production rate,
   exposure evaluator, accepted arc, event and allocator. The fixture is
   verification input. Choose its horizon and distance variation by
   independent analysis before inspecting AMPS values. Require numerical
   error below 0.2% with the approved strictly decreasing or exact-zero
   behavior. If a path is integrated analytically by the implementation,
   commission a separate force-driven path with step-dependent trajectory
   error and an independent dense trajectory/exposure reference. Preserve
   any new failures.
3. **Test physical geometry and all configured controls.** Run quiet and
   active D04 endpoints at prescribed 1 and 2 AU and confirm the rate ratio
   of four within the existing 0.1% tolerance. Exercise production Moon
   and Earth shadow services separately with analytically referenced
   finite-ray ingress and egress, zero dark hazard, and accepted
   optical-depth carryover through shadow. Run actual SPICE geometry cases
   separately with the existing hashed kernel set and named epochs. At one
   common physical horizon, compare the two force-aware movers, 8/16/32
   steps, 1/8 MPI ranks, and 1:1, 1:8 and 8:1 parent/product reference
   weights with fractional correction. Verify weighted nuclei, charge,
   mass, momentum and energy, bath transfers, and event ordering and
   ownership across cells and ranks. Use the already approved Na P02 limits
   in their named scope.
4. **Freeze the powered case matrix.** Use `I21/prepare_cohort.py` and the
   authenticated `m3_qualification_executor.py` entry path to generate and
   retain repeatable case receipts. Extend them and the existing I21 auditor
   to consume the trajectory-exposure cases. Size statistically independent
   parent cohorts with the approved planner. Score lifetime error ≤2%,
   inverse-square error ≤0.1%, and conversion budget closure after U14 and
   the applicable linked predecessors pass. [I21 execution contract](integration-tests/I21/EXECUTION_CONTRACT.md)

## 3. I22 — Commission atomic events, then close the molecular and surface gates

1. **Atomic pilot.** Start with ⁴He, then ²⁰Ne and ⁴⁰Ar, using the prepared
   isotope-specific ground-ion references. Keep the tabulated D04
   threshold, physical formation energy and source-state compatibility
   separate. Where the reviewed crosswalk supports it, pair the D04
   printed threshold with its quiet/active source mean to form the declared
   effective photon in the selected spectral convention. Use the
   physical mass/Q for the products. Test the actual event transaction
   with the existing photon bath/electron ledger and independently computed
   two-body recoil. A ground-ion reduced case needs a stated energy/state
   approximation and allowed error. The state reference alone supplies no
   branching fraction into that level. [Atomic family contract](../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
2. **Molecular reference profile.** Recommend para-H₂¹⁶O in
   `v = 0, J_KaKc = 0_00` as the first explicit energy/conservation
   reference. A separate ortho `1_01` sensitivity case can follow when
   both states and their ionic crosswalk are pinned. Freeze the actual OH
   and product reference levels, isotope masses and uncertainties, all on
   the P02 separated neutral-atom/electron energy reference. This choice
   supplies neither a lunar ortho/para ratio nor state-specific D04 rates.
   Keep O(³P), O(¹D) and O(¹S) distinct and check all populated-state
   exits. Use a prescribed zero-density collider bath for the initial
   reference and qualify spontaneous cascades. Require state- and
   density-dependent bounds for nonzero collisional quenching. The
   reported O(³P)+2H product assignment at 1216 Å provides a
   wavelength-specific cross-check. [Slanger and Black](https://doi.org/10.1063/1.444111)
3. **Missing source thresholds.** For molecular D04 channels without a
   printed numerical threshold, propose an explicit opt-in
   physical-reference-threshold plus source-mean effective-photon
   convention for controlled verification and sensitivity. Keep its
   reference-state identifier and both energy terms in every receipt.
   Admit each solar endpoint only when the source definition, state/Q,
   available energy and observable-specific error scope are compatible.
   Retain the D04 rates, means and seven/four-channel denominators.
   Preserve the problematic low-photon-energy OH tail in the archive.
   Wavelength-resolved use needs a documented interpretation of that
   tail. [Molecular rate and state inventory](M3_MOLECULAR_KINETICS_PROPOSAL.md)
4. **Charged breakup.** Pursue H₂O04 (OH⁺ + H + e⁻) first. The Bodi et al.
   threshold and breakdown curves anchor a near-threshold check. For the
   solar-conditioned closure, seek matching coincidence or evaluated
   differential data for the electron energy, the heavy-pair
   kinetic-energy release, product states and angular dependence. Use
   the same admission for H₂O05 (O⁺ + H₂ + e⁻) and H₂O06
   (H⁺ + OH + e⁻). A measured conditional distribution or a
   source-matched reduced law with independent error bounds can supply
   the missing partition. The numerical construction uses two Jacobi
   relative momenta. Subtract the selected internal energies and
   center-of-mass recoil, apply the sourced partition to electron versus
   heavy-pair motion, and allocate the remainder to ion–neutral relative
   motion with the corresponding reduced mass. Verify independent
   energy/momentum reconstruction, reversible frame transforms, and
   one-photon and one-electron exchanges in the actual transaction.
   Retain the approved neutral three-body law in its existing scope.
   Maintain all eleven source branches in audits and leave full I22
   SKIPPED while any required partition is unsupported.
   [Bodi et al.](https://doi.org/10.1039/C4SC00791C)
5. **U24 surface reference.** Start with the Laffon table’s 1 cm bare
   olivine sample and 1 cm gold control at 20 ± 1 K with 300 K background
   gas. Match the exposure sequence and compare the simulated
   sticking observable relative to the gold control. The olivine row
   reports `0.65 ± 0.07` relative to gold assigned unity. Retain the
   source's exposure-resolved values and uncertainty convention when
   freezing the comparison. Require separate evidence for absolute
   sticking or thermal accommodation. For water
   release, propose the published heterogeneous-site model with its
   documented material, mobility and coverage regime. Verify
   constant-temperature depletion independently before the coupled
   sticking/release calculation. For OH release, test the final
   Minissale ASW/olivine first-order pair together
   (`E/kB = 5698 K`, `ν = 3.76e15 s⁻¹`) in its source-defined setting.
   Complete incoming-OH sticking/accommodation and any required material
   mapping with separate evidence before U24 can pass in full.
   [Laffon source table](https://doi.org/10.6084/m9.figshare.13274840.v1)
   [Schörghofer water-release model](https://doi.org/10.3847/PSJ/ade5b2)
   [Minissale et al.](https://doi.org/10.1021/acsearthspacechem.1c00357)
6. **Full acceptance.** At both D04 endpoints, score I22 against
   independently computed source totals and channel probabilities with
   the entire frozen family. Retain absorption outside the reactive
   denominator. Apply the 2% ensemble rate limit, 0.2% numerical budget,
   95% linked branch consistency and approved branch precision separately.
   Retain the U15 99% criterion. Include molecular states, cascades, weights
   and exact represented nuclei/charge accounting in the independent audit.
   [I22 contract](integration-tests/I22/reference/acceptance.json)

## 4. I23 — Start with a bounded electron-impact helium pilot

1. Use an unpolarized stationary ⁴He target with a prescribed isotropic
   Maxwellian in the target rest frame. Recommend a first `Te = 10 eV`
   pilot at `ne = 0, 10⁶, 10⁷ m⁻³`, with the applicable D05 fit and
   physical ground-state threshold pinned independently. These explicit
   verification inputs test zero exposure and linear density scaling.
   Use the frozen full species and temperature plan for qualification.
2. In the authorized inherited-heavy-velocity model, evaluate the
   source-conditioned mean incident electron energy as
   `kB*Te*(3/2 + d ln k/d ln Te)` when the source/cross-section and
   unshifted-Maxwellian conditions of that identity apply. Bind the
   existing rate and moment services through the normal mover and
   chemistry transaction. Record one incoming electron, two outgoing
   electrons, the ion, the threshold and signed bath momenta and energies.
   With outgoing minus incoming signs the electron bath gains one
   electron and charge `−e` per event. Audit survival, represented
   budgets and chronology independently. A zero-density run retains all
   parents. A mismatched threshold, invalid plasma or unsupported drift
   must be reported explicitly.
3. Bound the omitted heavy recoil using the approved reservoir method,
   with `S = sqrt(2*me*Ein) + sqrt(4*me*Eout)` as an RMS impulse
   bound within its stated model. Report RMS speed bound `S/mi` and the
   associated absolute kinetic-energy discrepancy bound for the declared
   lab speed domain. Derive source/fit uncertainty for the conditional
   mean separately from the scalar rate error. Keep these approximation
   limits separate from the 0.3% numerical allocation. Reuse the service
   for Ne, Ar and Na only after each source-domain and observable
   crosswalk passes review. For Na, address the fit’s rounded 5.1 eV
   threshold versus the higher physical threshold through independent
   source evidence or a separately reviewed alternative provider. Preserve
   the D05 scalar rate verification.
4. Run the full qualified cases at ≤2% direct coefficient and ≤3%
   ensemble lifetime error. Require D06/R4.1 only when claiming measured
   plasma forcing. Continue to keep unsupported transport outside the
   admitted bath/target domain.
   [Reservoir contract](M3_RESERVOIR_PRODUCTS_CONTRACT.md)
   [I23 contract](integration-tests/I23/reference/acceptance.json)

## 5. I24 — Admit one defensible first transport channel

1. **First choice.** Time-box one source-review session to recover the
   original experimental citation and preparation, native relative-energy
   convention, uncertainty, prompt versus eventual H(1s) interpretation,
   and angular support for D119990. Keep its 12 native knots over
   4–25 keV/u. Admit proton + Na(3s) → Na⁺ + H(1s) only when the selected
   partial rate and its product states are mutually supported. With that
   crosswalk, independently check `Q = IE(H) − IE(Na)` using the named
   states. The source’s missing quantities keep this pilot gated.
   [IAEA CollisionDB D119990](https://db-amdis.org/collisiondb/datasets/119990/)
2. **Fallback order.** First seek the numerical state-to-state proton–Na
   data underlying Courbin-Gaussorgues et al. at 0.5–5 keV, with separate
   Na(3s) and Na(3p) targets, calculation uncertainty, and angular
   coverage or an independently justified bound. If that cannot be
   admitted, investigate a specific published proton–He energy/state
   pair in the Schöffler et al. coincidence measurements. Those
   measurements cover a much higher energy regime. Validate the
   energy convention and relativistic error for the selected case
   independently. Either fallback needs a case selection and scope
   review before replacing a frozen case. Keep the original ORNL
   helium and Morgan sodium total-capture source checks. A partial
   channel needs its own normalization and coverage proof before it can
   represent a total-capture model.
   [Courbin-Gaussorgues et al.](https://doi.org/10.1088/0953-4075/23/21/023)
   [Schöffler et al.](https://www.atom.uni-frankfurt.de/publications/files/Schoeffler2009pra.pdf)
3. **Production and oracle.** Start with a declared monoenergetic bath.
   Evaluate `n_projectile*σ(g)*g` from the actual projectile density
   and relative speed in the declared frame. Compare against an
   independent evaluation at source knots and in-range points.
   Verify the sampled collision-weighted projectile distribution before
   extending to a broad bath. Use conservative two-heavy-body kinematics
   with an evidenced differential law or a bounded forward approximation.
   For the latter, `2*p_final` bounds the omitted impulse over all angles.
   A smaller `2*p_final*sin(θmax/2)` needs supported angular coverage
   and treatment of any omitted tail. Check the resulting effect on the
   declared heavy-product speed and energy observables before enabling
   the model. Record the incoming prescribed ion and outgoing fast
   neutral separately. Use source-valid domains and retain unsupported
   channels OFF. Qualification retains ≤1% rate, ≤3% lifetime, 0.3%
   numerical bias, exact applicable budgets and the U17 minimum of
   100,000 probabilistic-branch events. [I24 contract](integration-tests/I24/reference/acceptance.json)

## 6. Execution, evidence and review sequence

| Phase | Deliverable | Exit condition |
| --- | --- | --- |
| One-hour development packages | I21 trajectory oracle and case plan first. In parallel, atomic I22 profile and I23 He provider/probe contracts. Use a separate bounded review package for the molecular/U24 and I24 source evidence. | Retained per-case results, explicit data/domain and approved-model limits, zero-byte nightly Moon diff after source changes. Unsupported cases retain their status. |
| Candidate linked pilots | Actual production events from admitted profiles for I21 and then I22/I23/I24 as each is ready. Commission common scorer/manifest support and the global extended adapter under the existing IDs. | Each case has an authenticated build and native receipt, independent numerical and budget evidence, and a measured cost. Default nightly coverage stays at its documented scope. |
| Cluster qualification | Freeze the complete statistical family, plan size, selections, resources, manifests and layout before submission. Include same-campaign I13–I20 and applicable U-series gates. | All original M3 acceptance conditions evaluated and met, with preservation of earlier scientific and core regression results. |

Use explicit one-hour packages (3,600 s including setup, build, execution,
plotting, reporting and cleanup), 8 MPI/OMP1 by default and an explicit
one-rank control where required. Enforce the 200 GB aggregate ceiling.
Allow reserve for preservation and cleanup. Archive observed wall time,
memory, disk and event counts, and record resource exits with the
limiting quantity and how to reduce the requested resource footprint.
The source-matched I21 pilot took 388.10 s for six launches and peaked at
42.93 GB sampled process-tree RAM. Size subsequent work from measured cost,
including MPI layout, parent count, branch rarity and build/setup overhead.
Use a separately budgeted cluster package for cases whose powered sample
cannot fit. Derive the CPU-hour request from measured rates and the full
frozen case plan before submission.

For reproducibility, run from a private source-matched checkout, remove its
`build/` before configuring, and use the normal Moon/SPICE configuration and
`make -j` build followed by `make -j test_Moon`. Require a zero-byte
`test_Moon.diff` and preserve applicable passing core regressions if shared
code changes. On the current host the build commands are as follows.

```bash
rm -rf build
./Config.pl -application=moon -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j
make -j test_Moon
test -f test_Moon.diff && test ! -s test_Moon.diff
```

Reuse the actual `MoonNativeReactionRegistration`,
`ReactionExecution`, force-aware mover and controller services. Define each
case in the existing U/I contract before code changes. Preserve disabled
legacy behavior and document new switches, units, frames, domains and
failures in code comments, READMEs and STATUS.

Freeze the plan and exact D04/D05 and optional D06 selection and hashes,
source catalog lock, source-state/error crosswalk, ordered SPICE identities,
case grid and chosen limits. The currently inspected shared chemistry
snapshot is `amps-chemistry-20261120-r5`, lock SHA-256
`7b588795ba3dcbd17491d5afbb9a417f62824970afa94b10645e8eda835799a2`.
Publish newly supported source material through the reusable chemical
database with immutable originals, versioned derived tables, tests and
licensing/provenance. Keep Moon-specific profile selections in the Moon
application. [Current source and data identities](M3_COMPLETION_PROPOSAL.md)

Use original-parent endpoint survival, the approved numerical allocations
and parent-keyed random streams at matched horizons across layouts. The primary I22
metric is total removal, with conditional branch and partial-rate diagnostics.
For I21/I23/I24 the primary metric is lifetime. Pre-size each active branch
for simultaneous 95% exact binomial intervals with relative half-width
≤10%. Preserve all seeds and zero-hazard outcomes. Record stdout/stderr,
`result.json`, `run_manifest.json`, input and configuration, build and
executable identity, data hashes, RNG and MPI layout, original-parent and
event/ledger records, comparison CSV/PNG/EPS and resource history under a
new `test_output/srcMoon/` campaign directory. Native completion,
preparation and scientific acceptance each keep their actual status.
[Approved scoring and planning details](M3_QUALIFICATION_ADDENDUM_20261118.md)

**Decisions to take first**

1. Approve the I21 varying-exposure cases and the rule for reviewing the
   constant cohort’s formal role with its retained failure visible.
2. Approve the I22 controlled atomic/para-water reference scope, the
   explicitly named effective-photon convention and the source-matched
   U24 reference strategy. Review the charged partition separately after
   a matching quantitative source is available.
3. Approve the stationary ⁴He I23 pilot at the proposed plasma inputs.
   Review the modeled observables and their source-supported error limits
   before physical activation or any extension to moving targets.
   *Resolved 2026-11-25 for the coupling scope:* the user approved the
   inputs and a narrow observable admission (heavy budgets and chronology).
   Moving targets and velocity-sensitive observables remain open.
   [Record](integration-tests/I23/HELIUM_ADMISSION_20261125.md)
4. Approve the I24 source-review order. Select the first transport case and
   its applicable energy/angle envelope after quantitative source review.
   *Resolved 2026-11-25:* the user selected total one-electron capture
   H⁺ + Na(3s) → Na⁺ + H(n) at 1–5 keV/u: the Morgan 1985 total rate, the
   Shingal 1986 state split, zero-angle exchange with a source-based recoil
   bound, and prompt H(n) with a radiative-decay energy ledger.
   The user approved the implementation details the same day, and the case
   is implemented with supporting ensemble and linked-PIC tests; formal I24
   remains SKIPPED.
   [Review and decision](integration-tests/I24/SOURCE_REVIEW_20261125.md) ·
   [Implementation contract](integration-tests/I24/TOTAL_CAPTURE_CONTRACT.md)

The immediate software work is case preparation, independent oracles and
runner/manifest integration. Physical activation and cluster qualification
depend on the applicable evidence, reviewed case scope and predecessor
results above.
