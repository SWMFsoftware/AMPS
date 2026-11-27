# I22 decision 2, parts B–D and charged breakup — proposals (2026-11-26)

**Proposal only, for the user's decision.** Nothing here is approved physics.
This work changed no production code, test, threshold, observation, raw data
or reference. Formal I22, U15 and U24 remain **SKIPPED / NOT QUALIFIED**.
Every value below comes from a local file or a source that was actually
opened, and each one cites that file or source. "Not found" means it was not
found; no value is filled in.

**Context.**
- Decision 2 of the [M3 follow-up](M3_COMPLETION_FOLLOWUP_20261124.md),
  section 3 and the decision list.
- Part A, the atomic pilot, is done ([STATUS 5.63](../STATUS.md)).
- Source: HEAD `01124778` (the atomic pilot commit).
- Roadmaps: the source copy (`712e1152…ee17`) and the data copy
  (`d10eca1c…49c4`) differ only in source-only amendments. None of those
  amendments touches U15, U24 or I22.

**Research evidence.**
- [`test_output/srcMoon/i22-decision2-research-20261126/`](../../test_output/srcMoon/i22-decision2-research-20261126/README.md)
  holds the downloaded originals with checksums.
- Earlier evidence used here:
  - [molecular rate inventory](M3_MOLECULAR_KINETICS_PROPOSAL.md);
  - [literature follow-up](M3_PHYSICS_LITERATURE_20261119.md);
  - [P02–P04 decisions](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md);
  - `test_output/chemical_database/m3-targeted-inputs-20261119T235459Z/molecular/report.md`,
    the threshold algebra.

## 0. Summary and recommended order

| Part | Proposal | New production work | What stays SKIPPED |
|---|---|---|---|
| **B** Molecular reference states | Para-H₂¹⁶O (0₀₀) and OH(X²Π₃/₂, v=0, J=3/2) parents on the P02 zero, with explicit product levels, O(³P)/O(¹D)/O(¹S) kept distinct, and a radiative cascade from one frozen A-value set. A **molecular ledger pilot** under I22: **OH family first**, which is complete with no charged breakup, then H₂O. | Molecular state profile; molecular photolysis provider; compiled H₂O/OH/O/H/H₂/ions at P02 masses; excited-O state tag with timed decay | Daughter kinematics, escape and brightness; ortho water; collisional quenching |
| **C** Missing thresholds | A named opt-in photon policy, "physical reference threshold + D04 mean". First comes an **evidence step (C-0)** that reconstructs the D04 rates and the implied PHIDRATES thresholds from the cross sections and the PHIDRATES solar flux. | A third `PhotonEnergyPolicy`, a molecular crosswalk, and an acquisition and processing script | Wavelength-resolved use and the interpretation of the OH O(¹S) tail |
| **D** U24 surface reference | Three supporting U24 sub-checks with no new ID: OH first-order desorption (Minissale), relative sticking (Laffon), and H₂O heterogeneous-site release (Schörghofer). The aggregate U24 stays SKIPPED. | OH/H₂O surface parameters; a coverage-dependent site inventory (for the Schörghofer check only) | Absolute sticking, incoming-OH accommodation, migration, cold traps and full U24 |
| **Charged breakup** H₂O04–06 | Staged: (0) audit the coverage of the D04 branch sources; (1) two-body H₂O⁺/OH⁺ ionization inside B; (2) a near-threshold H₂O04 controlled verification; (3) a solar-conditioned partition only once sources are in hand. While unsupported, H₂O04–06 events go to an explicit **untransported reservoir**. | Molecular wrapper for the existing `PhotonChargedBreakup` (test-only today) | Full H₂O-family I22 while any partition is unsupported |

**Recommended order:**
1. **C-0 evidence step.** It needs approval to acquire the PHIDRATES solar
   flux. It serves both C and the charged-breakup audit.
2. **The B + C pilot.** Run the OH family first, then H₂O with the reservoir
   for H₂O04–06.
3. **The D sub-checks**, in parallel. Start with OH first-order desorption,
   because its production kernel already exists.
4. **Charged-breakup stage 2.** Stage 3 waits for sources.

## 1. Common energy reference (B, C and charged breakup)

P02 sets the zero: separated neutral ground-isotope atoms H(1s ²S₁/₂) and
O(2p⁴ ³P₂) and a free electron at rest. From
[P02](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md):

- `m_ref(s) = Σ N_a m_ref(a) − z m_e`
- `E_model = |p|²/2m_ref + ε_formation + ε_excitation`
- `Q = Σε_in − Σε_out`, applying no excitation twice.

The formation energies are F(OH) = −Do, F(H₂O) = −(Dw+Do) and F(H₂) = −D₂.
For an ion, F(X⁺) = F(X) + IE(X).

**Primitive values found.** Uncertainties are as printed by each source;
"n.p." means not printed.

| Symbol | Quantity | Value | Source | Status |
|---|---|---|---|---|
| Dw | D₀(H₂O → H + OH), lowest | 41145.92 ± 0.12 cm⁻¹ (95 %) | Boyarkin et al. 2013, doi:10.1016/j.cplett.2013.03.007 (via NIST JRES 121.005, Burgess 2016, local lossless table) | Verified 2026-11-26 from the paper (user-supplied, `fba516bc…`). D₀ = E₀ + hν₁ + hν₂ + hν₃ᶜ is the continuum onset in total rovibrational energy. E₀ is the rotational energy of the initial level, so D₀ is measured from the 0₀₀ zero. The products are OH[X²Π₃/₂, e, J=3/2, v=0] + H(ground) |
| Do | D₀(OH → O + H) | 35584 ± 10 cm⁻¹ (Joens 2001); 35580 ± 15 cm⁻¹ (Sun et al. 2020); 35580 ± 3 cm⁻¹ (ATcT, as cited by Sun et al. 2020, ref. 46 = Ruscic 2015, doi:10.1021/acs.jpca.5b01346) | Joens 2001, doi:10.1021/jp011833u: a thermochemical cycle, D₀(H–OH) = AE₀(OH⁺/H₂O) − IE(OH), with 2001 inputs. Sun et al. 2020, doi:10.1063/1674-0068/cjcp2003026: HRTOF measurement. Both user-supplied | **Needs your choice.** Recommended: ATcT 35580 ± 3, best confirmed directly from ATcT; otherwise Sun 35580 ± 15. The 4 cm⁻¹ (0.5 meV) spread changes no threshold beyond its stated uncertainty |
| Io | IE(OH), X³Σ⁻ v⁺=0 | 104989 ± 2 cm⁻¹ = 13.01698 eV | Wiedmann et al. 1992 (ZEKE); Harper 2019 gives 13.01698(25) eV | Consistent |
| Iw | IE(H₂O), adiabatic | 101766.3(10) cm⁻¹ (Lauzin, Jacovella & Merkt 2015); 101766.8(1.2) and 101766(2) cm⁻¹ (Burgess 2016); 12.621 ± 0.002 eV (NIST WebBook evaluation, = 101795 cm⁻¹) | as listed | **The WebBook value differs from the spectroscopic values by ≈3.5 meV: needs your choice** |
| D₂ | D₀(H₂), para-H₂ X(v=0, N=0) | 36118.06962(37) cm⁻¹ | Liu et al. 2009, J. Chem. Phys. 130, 174306 (user-supplied PDF `d79e96e2…`; the abstract value and the para definition are on line 668 of the text extraction) | Verified 2026-11-26 |
| IH | IE(H) | 109678.77174307(10) cm⁻¹ | Burgess 2016 | — |
| IO | IE(O) → O⁺(⁴S₃/₂) | 109837.02(6) cm⁻¹ | NIST ASD / Burgess 2016 | — |
| ED, ES | O(¹D₂), O(¹S₀) above ³P₂ | 15867.862, 33792.583 cm⁻¹ (n.p.) | NIST ASD v5.12; same values in CHIANTI 11.0.1 `o_1.elvlc` (`230bbdcd…`) | Consistent |
| — | Ortho 1₀₁ − para 0₀₀ | 23.79436122(25) cm⁻¹ | Tóbiás et al. 2020, doi:10.1038/s41467-020-15430-6 | Only for a later ortho case |

**Derived reference thresholds** (eV, CODATA 2018 `hc`). The ± column adds
the printed uncertainties in quadrature, which assumes independence;
covariance was not evaluated.

| Channel | Algebra | Threshold [eV] | Check |
|---|---|---|---|
| H₂O01 H + OH | Dw | 5.10144 | — |
| H₂O02 H₂ + O(¹D) | Dw+Do−D₂+ED | 7.00259 ± 0.0012 | — |
| H₂O03 O(³P) + 2H | Dw+Do | 9.5133 ± 0.0012 | — |
| H₂O04 OH⁺ + H + e | Dw+Io | 18.11842 ± 0.00025 | **Bodi 2014 measured 18.1182 ± 0.0009 eV; agrees** |
| H₂O05 O⁺ + H₂ + e | Dw+Do−D₂+IO | 18.65328 ± 0.0012 | Kronebusch & Berkowitz 1976, Table 7: 18.65 eV |
| (O⁺ + 2H + e, atomization) | Dw+Do+IO | 23.13135 ± 0.0012 | Kronebusch & Berkowitz 1976, Table 7: 23.126 eV |
| H₂O06 H⁺ + OH + e | Dw+IH | 18.69988 | Sage 2010 D₀(H⁺–OH) route gives 18.697–18.700 eV |
| H₂O07 H₂O⁺ + e | Iw | 12.61741 ± 0.00012 (Lauzin) | — |
| OH01 O(³P) + H | Do | 4.4119 ± 0.0012 | — |
| OH02 O(¹S) + H | Do+ES | 8.6016 ± 0.0012 | — |
| OH03 OH⁺ + e | Io | 13.01698 ± 0.00025 | — |
| OH04 O(¹D) + H | Do+ED | 6.3792 ± 0.0012 | — |

ATcT (Ruscic) was the first-named reference family in P02. Every request to
`atct.anl.gov` returned HTTP 403, so no ATcT value is used here.

## 2. Part B — molecular reference states

### 2.1 Proposed reference states

| Species | Proposed reference level | Basis |
|---|---|---|
| H₂O parent | para-H₂¹⁶O, X̃ ¹A₁, v=0, J_KaKc=0₀₀ (the absolute lowest level) | Follow-up recommendation. Ortho 1₀₁ is a later sensitivity case |
| OH parent and product | ¹⁶OH X²Π₃/₂, v=0, J=3/2 (lowest) | The usual zero for D₀(OH) |
| H₂ product (H₂O02) | para-H₂ v=0, J=0 | Lowest level. Proton nuclear-spin symmetry makes para water correlate with para H₂; **please confirm this physics reasoning** |
| O | ³P₂ / ¹D₂ / ¹S₀, kept as distinct states | D04 labels: H₂O02 and OH04 O(¹D); OH02 O(¹S); OH01 O(³P) (rate HTML) |
| O in H₂O03 | **Unlabelled in D04** ("O + H + H") | Option: declare O(³P) as a reduced choice, with Slanger & Black 1982 at 1216 Å as support. Record O(¹D)+2H, which Chang et al. 2021 observe at 92–109 nm, as an eligible but unpopulated state. This follows the approved Ne⁺/Ar⁺ handling |
| Ions | H₂O⁺ X̃ ²B₁ (0₀₀); OH⁺ X³Σ⁻ v=0, N=0; O⁺ ⁴S₃/₂; H⁺ | Ground ion, as in the atomic pilot. Ã/B̃ ion populations are recorded as eligible, not as yields |
| Masses | P02 reference masses from the NIST isotope masses already in the atomic profile (`1H:ground`, `16O:ground`) | `input/species.input` uses integer AMU (`18*_AMU_`); it must not be used for the ledger |

What D04 does and does not supply:
- The D04 rates are for the laboratory cross-section sample, a thermal
  population. They are not state-resolved, so choosing para 0₀₀ fixes the
  energy-ledger zero, not the rate.
- The scalar mean excess E does not fix any daughter distribution. This is
  recorded in the kinetics proposal.
- The current atomic-family `states.csv` already lists the molecular and
  excited-O rows as `unavailable`. These include `H2-16O:unresolved`,
  `16OH:unresolved` and `16O:1D-unresolved`.
- The profile loader refuses a molecule on the atomic zero
  (`molecule-cannot-use-atomic-formation-zero`). A molecular profile
  therefore needs actual `energy`-role sources, following
  `ATOMIC_FAMILY_CONTRACT.md`: authenticated originals, printed wavenumbers
  with uncertainties, exact SI conversion, and a native registry probe.

### 2.2 Excited oxygen

Both A-value sets are local or opened. They disagree:

| Set | τ(¹D₂) | τ(¹S₀) | Note |
|---|---|---|---|
| CHIANTI 11.0.1 `o_1.wgfa` (`f24efc2b…`; Froese Fischer & Tachiev 2004) | **116.22 s** (ΣA = 8.6046e-3 s⁻¹) | **0.74695 s** (ΣA = 1.33877 s⁻¹; 5577 Å fraction 0.9412) | Recomputed here from the native file. The Nov-19 review found its green/297 nm photon ratio (16.07) close to the Mars observations (16.5 ± 0.4; 15.8 ± 2.4) |
| NIST ASD v5.12 | **133.8 s** | 0.7487 s | No uncertainties printed |

The production service `srcMoon/species/OxygenCascade.{h,cpp}` already
implements a five-level, nine-link CHIANTI network, but only tests call it.

Proposal:
- Use one frozen A-value set; the recommendation is CHIANTI 11.0.1.
- Treat radiative decay as a competing timed event that keeps the atom and
  its trajectory, updates its state, and books the photon energy and
  momentum.
- Use a zero-density collider bath, so there is no quenching in the pilot.
  Any nonzero quenching needs separate state- and density-dependent bounds.

### 2.3 Proposed test: I22 molecular ledger pilot (existing ID, supporting only)

- **Purpose.** Run complete D04 families as real particle events through the
  normal cohort owner, as the atomic pilot does. Verify chronology, cascade
  timing and an exactly closed nucleus, charge, mass, energy and momentum
  ledger against independent oracles.
- **Production capability needed.**
  - A molecular photolysis provider next to `AtomicPhotoionization`, binding
    the complete D04 family of the parent.
  - Kinematics:
    - neutral two-body: isotropic in the parent rest frame, with `K_avail`
      from the Part C policy;
    - two-body ionization: the P02/P03 photoelectron treatment;
    - H₂O03: the approved P04 three-neutral sampler
      (`ReducedKinematics::ThreeNeutralMicrocanonical`);
    - excited-O decay events.
- **Entry points.**
  - `Moon::Species::PhysicalReactionEvents`, which needs a new molecular
    target;
  - `Moon::Chemistry::SodiumCohort::Prepared`, the generic cohort owner, and
    the native adapter;
  - `ReducedKinematics`, `OxygenCascade` and the force-cut-surface movers.
- **Configuration.**
  - A zero-acceleration compiled build with the parent and all products at
    P02 masses.
  - Frozen LSO and SPICE Sun as in Part A.
  - Quiet and active endpoints; both movers; 8/16/32 steps; 1 and 8 ranks;
    seed 731993.
  - The substep cap equals the coarsest step (Part A revision b).
- **Domain constraint (new).** Under the mean-energy reduced model, H
  daughters reach about 4.2×10⁴ m/s (OH02: 9.4 eV goes to H). That
  limits the horizon to roughly 7×10² s from the current start
  point if every product must stay inside the domain with the Part A
  1.1×10⁷ m margin (5.21×10⁷ − 10⁷ − 1.1×10⁷ m at 4.2×10⁴ m/s ≈ 730 s).
  Two alternatives:
  - a domain-centred, off-shadow start point, where |y|,|z| set the limit at
    about 1,330 s;
  - a new domain-exit ledger category, which is a production and verifier
    change.
  With the centred start, a 1,024 s horizon gives an OH-family removal
  probability of about 2.0 % (quiet) and 3.4 % (active). About 5,000 and
  3,000 parents then give about 100 events.
- **Oracles.**
  - (1) An exact event-ledger reconstruction counting H and O nuclei, with
    multiplicity 2 for H in H₂O03, plus charge, mass, energy and momentum.
  - (2) The straight-line integrated-exposure event-time oracle.
  - (3) Independent multinomial channel frequencies against `k_i/Σk`.
  - (4) The excited-O populations against a matrix exponential of the frozen
    radiative-rate matrix (the Bateman solution).
  - (5) Original-parent identity between 1 and 8 ranks.
  - (6) A verifier-side target table independent of the C++ code.
- **Metrics and acceptance.** The existing limits, unchanged:
  - event-time refinement strictly decreasing and ≤ 0.002;
  - exact ledger closure within the declared arithmetic bounds;
  - 1/8-rank identity;
  - channel frequencies within the existing U15 99 % familywise intervals;
  - cascade populations within binomial 99 % intervals of the
    matrix-exponential prediction. This last tolerance is proposed and
    needs your approval before any run.
- **Admitted observables.**
  - Admitted: heavy nucleus, charge and mass budgets; event chronology;
    excited-O state chronology.
  - **Not admitted:** daughter speed or energy distributions, escape,
    brightness, and state-resolved yields.
- **Expected failures.**
  - A product leaving the domain (ERROR).
  - An unmatched or missing state in the profile (SKIPPED).
  - Ledger non-closure, a frequency outside its interval, or non-decreasing
    refinement (FAIL).
- **Data.** D04 (`photo_rates_1au.csv` `3aa54915…`, `reaction_channels.csv`
  `2d9fa1c8…`), the molecular profile (new), and CHIANTI `o_1.*`.
- **Predecessors.** Recorded, not blocking, under the M3 development
  exception: U15, U24 and same-campaign I13–I20.
- **Status.** A supporting PASS at most. Formal I22 stays SKIPPED.

### 2.4 Decisions for Part B

- **B1.** Parent references: para-H₂O 0₀₀ and OH X²Π₃/₂ J=3/2 v=0
  (recommended), with ortho later.
- **B2.** Primitive values:
  - Dw: Boyarkin 41145.92(12) cm⁻¹, from 0₀₀ to OH(X²Π₃/₂, e, J=3/2,
    v=0) + H (verified).
  - Do: ATcT 35580(3) (recommended), Sun 35580(15) or Joens 35584(10) cm⁻¹.
  - Iw: Lauzin 2015 (recommended) or the WebBook value.
  - Para D₂: Liu et al. 2009, 36118.06962(37) cm⁻¹ (acquired 2026-11-26).
- **B3.** Product references as in 2.1. Also the H₂O03 oxygen state:
  - declare O(³P) with O(¹D)+2H recorded as eligible (recommended), or
  - leave it unresolved, which leaves the H₂O family incomplete and so not
    runnable.
- **B4.** A-value set: CHIANTI 11.0.1 (recommended) or NIST ASD v5.12.
- **B5.** A zero-density collider bath for the pilot.
- **B6.** The admitted observables and the cascade-population tolerance
  in 2.3.
- **B7.** Domain handling: short horizons with a centred start
  (recommended, in the spirit of Part A), or a new domain-exit category.
- **B8.** Order: the OH family first, then H₂O.
- **B9.** P02 reference masses for the compiled molecular species in the
  pilot build.

## 3. Part C — missing source thresholds

### 3.1 The problem

`thresholds.csv` (`a3335b3b…`) has rows only for He, Na, Ne and Ar. The
H₂O and OH_exp rows are absent. The raw PHIDRATES files print no numerical
molecular threshold either:
- the H₂O BrnOut blocks carry only `Lambda = …` ranges;
- the OH_exp total block names "Threshold: Huber & Herzberg 1979" with no
  number.

The D04 mean excess E is defined relative to the threshold PHIDRATES itself
used (Huebner & Mukherjee 2015). The production code refuses any channel
without a printed threshold (`PhysicalReactionEvents.cpp:255`). Inferring a
threshold from the first nonzero grid point is not allowed.

### 3.2 Proposed convention

Add a named, opt-in `PhotonEnergyPolicy`, `PhysicalThresholdSourceMean`:
- E_γ,eff = T_phys + E_D04, where T_phys comes from section 1.
- K_avail = E_D04 − photon recoil, which follows from E_γ,eff − T_phys.
- Every receipt keeps T_phys with its source identity, E_D04 with its row
  identity, and the fact that the PHIDRATES threshold T_src is **unknown**.

What it claims and does not claim:
- The D04 rates are unchanged.
- The error in K_avail equals T_src − T_phys, which is unknown until C-0.
- It is not admitted for wavelength-resolved use.

### 3.3 Evidence step C-0 (proposed before adopting the policy)

1. **Acquire** the PHIDRATES binned solar flux for the quiet and active
   endpoints from the public PHIDRATES site. Store it in a new versioned
   package written by a committed script; raw data stay immutable.
2. **Reconstruct** each channel rate, k_i = ∫σ_i φ dλ, from the local
   `cross_sections.csv` (`e900eb55…`). Compare it with
   `photo_rates_1au.csv`. The agreement tolerance must be declared before
   the run; I propose 1 %, report-only unless you approve it as a gate.
3. **Compute** the rate-weighted mean photon energy ⟨E_γ⟩_i. Then the
   implied PHIDRATES threshold is T_src,i = ⟨E_γ⟩_i − E_D04,i, and the
   policy error is ΔT_i = T_src,i − T_phys,i.
4. **Quantify, using the same integrals:**
   - the share of each H₂O04–06 rate lying outside 304–690 Å, the only range
     with a cited branch reference (Kronebusch & Berkowitz 1976);
   - the rate contributed by the OH_exp O(¹S) tail, which runs to 6132 nm,
     where the photon is about 0.2 eV and below the O(¹S) channel threshold
     (8.60 eV).

### 3.4 Electron speed and the classical-energy budget

These figures assume the electron takes the whole D04 mean.

| Channel | Mean (quiet/active) [eV] | v_e [m/s] | ¾(v/c)² |
|---|---|---|---|
| H₂O07, OH03 | 12.4/15.2; 19.4/23.6 | ≤ 2.88×10⁶ | ≤ 6.9×10⁻⁵ |
| H₂O04, H₂O06 | 18.6/23.2; 25.0/30.5 | ≤ 3.28×10⁶ | ≤ 9.0×10⁻⁵ |
| **H₂O05** | **36.5/39.8** | **3.58/3.74×10⁶** | **1.07/1.17×10⁻⁴** |

H₂O05 exceeds the I22 cap of 3.4×10⁶ m/s and the P02 relative limit of
1e-4. The options are:
- exact relativistic electron bookkeeping in the bath;
- the reservoir for H₂O05 (section 5.4);
- a new molecular cap with its own budget. This option would need explicit
  approval; it is not a relaxation of an existing gate.

### 3.5 Decisions for Part C

- **C1.** Adopt the opt-in policy `PhysicalThresholdSourceMean` with the
  receipts in 3.2.
- **C2.** Approve evidence step C-0, including acquiring the PHIDRATES solar
  flux and the report-only 1 % rate-reconstruction tolerance or another one.
- **C3.** The error scope: budgets and chronology do not depend on T, and
  kinematic observables stay not admitted.
- **C4.** How to handle the H₂O05 electron speed.
- **C5.** Record the OH_exp O(¹S) tail as an archived source feature, without
  clipping or repairing it.

## 4. Part D — U24 surface reference

### 4.1 Current state

U24 is "H2O/OH surface-migration and photochemistry kernel verification"
(roadmap R3.2, R4.7 and §8.1; same text in both roadmaps).
- **Contract.** `U24_water_chemistry/reference/acceptance.json` says
  acceptance is "Frozen only after physical parameters and production
  channels are documented". So **U24 has no numerical criterion yet**.
- **Status.** SKIPPED: "H2O and OH are not production species and no
  migration/chemistry path exists".
- **Only child check.** The production first-order Arrhenius
  `ResidenceTimeS`/`DesorptionProbability` against the Laferriere law,
  9/9 PASS.
- **Production surface physics.** It covers `Species {Sodium, Argon}` only.
  Sticking is a temperature table, with no material, coverage, exposure or
  angle dependence. There is no heterogeneous-site law. R4.7's
  `MoonWater.h`, `MoonPhotochemistry.h` and `MoonColdTrap.h` do not exist.

### 4.2 Sources

| Source | Local | Content used |
|---|---|---|
| Laffon et al. 2021, Figshare 13274840.v1 (CC BY 4.0) | Yes, `904badaf…`, CHEMICAL_DATABASE r14 | "H2O at 20K": Gold 1.0 cm, 1 at every exposure (S = 1 ± 0.07). Olivine 1.0 cm: 0.56, 0.52, 0.67, 0.81, 0.71 at 0.1–1 L (S = 0.65 ± 0.07). The measurement is XPS O1s uptake with isotropic 300 K background dosing; S is the mean over the five exposures |
| Minissale et al. 2022 (final PDF; "Taverne" open access, not CC) | Yes, `0017f390…`; evaluation permission **unknown** | OH on compact ASW and olivine, sub-monolayer: recommended ν = 3.76×10¹⁵ s⁻¹ and E = 5698 K (TST prefactor), T_peak 140 K, stated precision 10 %. KIDA lists the same pair |
| Schörghofer 2023 (PSJ 4:164) and 2025 (PSJ 6:164), CC BY 4.0 | Local HTML is a captcha page. Public pages were read, and the model code `desorptionrates.f90` (v1.2.3) is local | Site-energy distribution E_p ≈ 0.65 eV, W ≈ 0.22 eV. Prefactor 1e16 s⁻¹ (2023) or 2.2×10¹⁷/√T (2025 code). Two endmembers: efficient diffusion and no lateral diffusion. The 2025 Eq. 4 and Eq. 9 were images that could not be read |
| Jones et al. 2020, Zenodo 3822365 (CC BY 4.0) | Yes, 14 files | Apollo 14163/10084 TPD and site distributions. **The energy unit is not declared** in the native headers |

### 4.3 Proposed supporting sub-checks (no new ID; aggregate stays SKIPPED)

- **U24-a: OH first-order desorption (Minissale pair).** This one is
  cheapest because the production kernel exists. It needs only OH surface
  parameters.
  - Oracle: closed-form constant-T depletion,
    N(t) = N₀·exp(−ν t e^(−E/kT)).
  - The source is weak for OH. The local Minissale review states that
    OH's "experimental evaluation of the desorption parameters [is] almost
    impossible" because OH reacts during the heating ramp. The 4600 K value
    is "reported but not detailed" (its ref. 235). The recommended pair is
    a transition-state-theory transformation of it. No OH heating rate
    exists, so a peak-temperature check would be circular and is dropped.
  - Acceptance: closed-form agreement at the existing child tolerance of
    1e-12 relative. This verifies the kernel and the transcription only.
  - Scope: the source-defined substrate and sub-monolayer regime only, not
    lunar regolith.
- **U24-b: Laffon relative sticking.**
  - Production sticking is an input. Reproducing 0.65 therefore verifies
    only transcription and wiring: the impact-counting and exposure
    operator. It is not physics.
  - The proposed check drives an isotropic 300 K Maxwellian flux through a
    thin adapter. It compares impacts per area against the Hertz–Knudsen
    oracle Φ = P/√(2πmkT), with a binomial 99 % interval. It also checks
    that the configured olivine/gold ratio equals 0.65.
  - The exposure-resolved values (0.52–0.81) cannot be represented without
    a coverage-dependent law, so they are recorded, not scored.
- **U24-c: H₂O heterogeneous-site release (Schörghofer).**
  - This needs a new production kernel: a coverage-dependent site
    inventory.
  - Oracles, sub-monolayer only (v ≤ 1):
    - efficient diffusion:
      v(t) = [v₀^(−γ) + γ a t]^(−1/γ), with γ = W/kT and
      a = ν e^(−E_p/kT)/(1+γ). This is the analytic solution of
      dθ/dt = −S with Schörghofer 2025 Eq. 4,
      S = θν e^(−E_p/kT)(θ/θₘ)^(W/kT)/(1 + W/kT) for θ ≤ θₘ, now verified
      from the paper. Eq. 9 (2025) gives ν = 2.2×10¹⁷/√T s⁻¹ (T in K),
      replacing the 2023 constant 10¹⁶ s⁻¹; E_p ≈ 0.65 eV, W ≈ 0.22 eV,
      θₘ = 10¹⁹ m⁻²;
    - no diffusion: θ(t) = ∫ n₀(E) exp(−ν t e^(−E/kT)) dE by independent
      quadrature.
  - The model code's multilayer branch computes S2 but returns S1. This
    should be reported, not used.
- **U24-d: incoming-OH sticking and accommodation.** No evidence was found,
  so it stays SKIPPED.

Full U24 also needs the R4.7 scope: H₂O/OH as production surface species,
hopping, ballistic migration, cold-trap inventory convergence and the global
H/O budget. All of that remains open.

### 4.4 Decisions for Part D

- **D1.** Approve U24-a, U24-b and U24-c as supporting sub-checks under U24
  (no new ID), with the tolerances above set before any run.
- **D2.** Each sub-check is scoped to its source material and regime only.
  No lunar-regolith mapping is implied.
- **D3.** Laffon's "± 0.07": absolute standard deviation, as printed beside
  S (proposed), or relative. The SI text quotes ±7 % for flat surfaces.
- **D4.** Schörghofer version: the 2025 code, ν = 2.2×10¹⁷/√T
  (machine-readable), or 2023, ν = 1e16. Also the endmember: start with
  efficient diffusion and add no-diffusion as a second case.
- **D5.** Licensing: may the Minissale PDF (Taverne) and the Schörghofer
  code ("evaluation permission unknown") be used for verification?
- **D6.** New production physics: OH/H₂O surface parameters for U24-a and b,
  and the site-inventory kernel for U24-c.
- **D7.** Whether to pursue the Jones 2020 Apollo TPD as a later
  real-specimen benchmark. Its undeclared energy unit must be resolved from
  the source, not assumed.

## 5. Charged breakup (H₂O04 OH⁺+H+e, H₂O05 O⁺+H₂+e, H₂O06 H⁺+OH+e)

### 5.1 What D04 says

| Channel | Rate quiet/active [s⁻¹] | Mean quiet/active [eV] |
|---|---|---|
| H₂O04 | 5.54e-8 / 1.51e-7 | 18.6 / 23.2 |
| H₂O05 | 5.85e-9 / 2.21e-8 | 36.5 / 39.8 |
| H₂O06 | 1.31e-8 / 4.07e-8 | 25.0 / 30.5 |

All three branch references are Kronebusch & Berkowitz 1976, "Lambda =
304.-690. A", that is about 18.0–40.8 eV (`brnout_001.raw`, `b1e9efbf…`).
The branch partials are nonzero down to 1 Å, beyond any cited branch range.

H₂O05's mean (36.5 eV) plus its threshold (> 18 eV) puts the rate-weighted
photon energy above 40.8 eV. So part of the rate comes from photon energies
with no cited branch source; C-0 quantifies how much.

The D₂O measurement of Hult Roos, Eland et al. 2018 concludes that O⁺ and D⁺
from the C̃ state come from **full atomization**, not O⁺ + D₂. The D04
branch source itself agrees. The Kronebusch & Berkowitz 1976 abstract says
that ionization from the inner (2s-like) valence state of H₂O "leads to …
O⁺ fragments … with probable shattering of the molecular ion into atomic
and ionic fragments". Their Table 7 lists both H₂ + O⁺ (18.65 eV) and
2H + O⁺ (23.126 eV) energetics. The D04 product label "O⁺ + H₂ + e"
therefore conflicts with its own source at those energies. Tan et al. 1978
Table 4 gives an O⁺ appearance potential of ≈ 23 eV, which matches the
atomization threshold (23.13 eV), not O⁺ + H₂ (18.65 eV). This must be
recorded as a source conflict (decision E5).

### 5.2 Evidence found

| Source | Numerical content | Form and access |
|---|---|---|
| Bodi et al. 2014 (local, CC BY-NC 3.0) | OH⁺ onset E₀ = 18.1182 ± 0.0009 eV (iPEPICO, < 2 meV resolution). Threshold in the B̃ ²B₂ band, leading to OH⁺(X³Σ⁻) | Breakdown curves plotted only (17.90–18.12 eV). No KER, angles, or H⁺/O⁺ onsets |
| Hult Roos, Eland et al. 2018 (local accepted manuscript), **D₂O** | B̃-state yields: D₂O⁺ 10 ± 2 %, OD⁺ 70 ± 3 %, D⁺ 19 ± 2 %, O⁺ 0.15 ± 0.05 %. About ⅓ of the excess energy goes to translation (TOF widths, "less precise"). C̃ state: OD⁺ 16 %, O⁺ 24 %, D⁺ 60 %; CM KER lower limits of 5, 3 and 4 eV | Spectra are figures; the yields are in text. Isotope substitution |
| Gallagher et al. 1988 JPCRD 17, 9 (local scan) | X/A/B/C state partials and β; fragment-ion cross sections (Tan 1978; Cairns 1971); ±15 % for state partials, ±10 % for fragments | **Plotted only.** The tables are in JILA Data Center Report 32, not found |
| Sage et al. 2010 (local) | H₂O⁺(B̃) → OH⁺(X): isotropic recoil; D₀(H–OH⁺) and D₀(H⁺–OH) | Ions prepared by photoexcitation, not by photoionization |
| Harper 2019, NIST WebBook | IE(OH), OH⁺ yield 12.6–15 eV; appearance energies OH⁺ 18.115 ± 0.008 eV, H⁺ 18.7 ± 0.05 eV (electron impact), O⁺ 19.0 eV | — |
| Tan et al. 1978 (user-supplied scan, `3e6d839d…`) | Table 2: photoabsorption and fragment-ion (H₂O⁺, OH⁺, H⁺, O⁺) oscillator strengths against energy, 6 eV upward. Table 3: partial oscillator strengths for the 1b₁, 3a₁, 1b₂ and 2a₁ orbitals, 16.5–60.5 eV. **Table 4 (H₂O): appearance potentials and fragment-ion kinetic energies.** OH⁺ 18.1 eV, thermal. O⁺ ≈ 23 eV, 0.10 ± 0.02 eV. H⁺ 18.7 eV, 0.14 ± 0.02 / 0.18 ± 0.03 / 0.34 ± 0.04 eV at 25 / 30 / 50 eV. Fig. 8: H⁺ attributed to the 1b₂ and 2a₁ states. Stated systematic errors about ±5 % | **Tabulated.** It is a scan with a garbled OCR layer, so the values need reviewed manual transcription (decision E3) |
| Kronebusch & Berkowitz 1976 (user-supplied scan, `31f79c84…`) | Table 7: H₂O fragmentation ratios at nine line-source energies, 16.76–40.80 eV, with energetics. Table 6: 21.2 eV intensities against earlier work | Tabulated; scan, same transcription rule |
| Eland 1975 (user-supplied scan, `6d37fde3…`) | **D₂O⁺** breakdown (D₂O⁺, OD⁺, D⁺) against ion internal energy, and KER. As summarized by Kronebusch & Berkowitz: OD⁺ KER about 30 % of the available energy | Figures; D₂O only (decision E4) |
| Powis & Reynolds 1991 (user-supplied, `35b93456…`) | **D₂O⁺** B̃ ²B₂: He I PEPICO, D⁺/OD⁺ branching and dissociation dynamics | D₂O only (decision E4) |
| Banna et al. 1986 (user-supplied, `1f4c823d…`) | Tables I–IV: β, relative and absolute (Mb) partial cross sections per water orbital, 30–140 eV | Tabulated, text layer |
| Truesdale et al. 1982 (user-supplied, `c5a88dde…`) | Tables I–II: photoelectron asymmetry parameters β for the water ionic states | Tabulated, text layer |
| Norwood, Ali & Ng 1991 (user-supplied, `05d0d978…`) | PEPICO breakdown diagrams for energy-selected **H₂O⁺** and D₂O⁺ (B̃ ²B₂), 625–700 Å (internal energy 129–166 kcal/mol). Table I gives relative ion intensities. The authors state the H₂O and D₂O diagrams are "essentially identical". 3–5 % stable H₂O⁺ on a ~10 µs scale | Tabulated (Table I) plus figures. Direct support for E4 for **B-state branching only**; it gives no KER distributions |
| Not acquired, low priority: Dutuit et al. 1985 (mainly neutral fragment fluorescence); Suárez et al. 2018 (theory) | — | — |

**Unsupported at solar-conditioned energies.** None of the following was
found in a form that is both tabulated and accessible:
- the photoelectron versus ion-level split;
- the H₂O state-to-fragment branching versus photon energy;
- H₂O KER distributions;
- product internal states;
- ion–neutral recoil anisotropy;
- the co-fragment identity at C̃-state energies.

### 5.3 Staged proposal

- **Stage 0 (audit).** Part of C-0. Measure the branch-source coverage, and
  record the O⁺/H⁺ co-fragment conflict and the 304–690 Å limit.
- **Stage 1 (inside the Part B pilot).** H₂O07 H₂O⁺ + e and OH03 OH⁺ + e as
  two-body ionization to the ground ion, like the atomic pilot.
  - The ion recoil from the electron is small (H₂O⁺ ≈ 64–70 m/s,
    OH⁺ ≈ 84–93 m/s), so heavy-ion transport does not depend on the
    unresolved ion-state partition.
  - Ã/B̃ ion populations are recorded as eligible, not as yields.
- **Stage 2 (controlled verification).** H₂O04 at a monochromatic photon
  between 18.1182 eV and the H⁺ + OH onset (≈ 18.70 eV), to OH⁺(X) + H + e.
  - It uses the existing `PhotonChargedBreakup` two-Jacobi construction with
    a **declared numerical partition**, labelled as numerical only.
  - Oracle: independent reconstruction of energy and momentum, reversible
    frame transforms, and one-photon/one-electron exchange in the actual
    transaction.
- **Stage 3 (solar-conditioned).** Sample the populated H₂O⁺ state from
  state partials, the photoelectron energy as photon energy minus the state
  energy, and the KER from measured distributions. This needs Tan 1978,
  Eland 1975 and Powis & Reynolds 1991 (a D₂O⁺ study), and a policy on digitization. Then
  H₂O05 and H₂O06 follow.
- **Interim, inside any H₂O-family run.** Execute H₂O04–06 as events and keep
  their hazards; never delete one and rescale the others.
  - Option (a), recommended: book their products to an explicit
    **untransported reservoir**. The reservoir keeps nuclei, charge, mass and
    energy; there are no product particles.
  - Option (b): the "inherit-parent-velocity, unresolved-electron reservoir"
    reduced model from the Nov-19 review. Its ions are transported, but its
    neutral fragments would carry the parent velocity, which is unphysical
    for H; e.g. a 0.5 eV KER gives H about 9.5 km/s.

### 5.4 Decisions for charged breakup

- **E1.** Approve the stages and the interim reservoir: option (a)
  (recommended) or (b).
- **E2.** May I obtain the paywalled papers through your institutional
  access (Tan 1978, Kronebusch & Berkowitz 1976, Eland 1975, Powis &
  Reynolds 1991, Banna 1986, Truesdale 1982)?
- **E3.** Digitizing plotted curves (Bodi, Gallagher, Hult Roos):
  - not allowed, so tabulated data only; or
  - allowed with a stated digitization-uncertainty budget and an
    independent second digitization.
- **E4.** Whether D₂O-to-H₂O transfer of the state branching and KER is
  acceptable, and with what bound. Until then it stays unsupported.
- **E5.** Record the D04 "O⁺ + H₂" label against the measured C̃-state
  atomization as a source conflict, without relabelling D04.
- **E6.** Stage 2's declared numerical partition. I would propose a value
  only after E3/E4, so none is assumed here.

## 6. Sources requested for acquisition (verified against Crossref, 2026-11-26)

The user supplied all thirteen papers below on 2026-11-26, in two batches. They are preserved
with SHA-256 in
`test_output/srcMoon/i22-decision2-research-20261126/user-supplied-20261126/`.
Later the same day the user also supplied Norwood 1991, Boyarkin 2013,
Schörghofer 2023 and 2025, Joens 2001 and Sun 2020, preserved in the same
place. With those, no further paper is required for these proposals. Still
needed: the PHIDRATES solar flux (data, decision C2). Optionally, a direct
ATcT confirmation of D₀(OH) = 35580 ± 3 cm⁻¹.

References for the second batch:
- Norwood, K., Ali, A., & Ng, C. Y. (1991). A photoelectron–photoion
  coincidence study of H₂O, D₂O, and (H₂O)₂. *J. Chem. Phys.*, 95(11),
  8029–8037. doi:10.1063/1.461334. It is the only identified H₂O
  state-selected source, so it reduces reliance on D₂O transfer.
- Boyarkin, O. V., Koshelev, M. A., Aseev, O., Maksyutenko, P., Rizzo,
  T. R., Zobov, N. F., Lodi, L., Tennyson, J., & Polyansky, O. L. (2013).
  Accurate bond dissociation energy of water determined by triple-resonance
  vibrational spectroscopy and ab initio calculations. *Chem. Phys. Lett.*,
  568–569, 14–20. doi:10.1016/j.cplett.2013.03.007. It fixes the
  initial and final levels of Dw.
- Schörghofer, N. (2023). Adsorption kinetics of water and argon on lunar
  grains. *Planet. Sci. J.*, 4(9), 164. doi:10.3847/PSJ/acf19b.
- Schörghofer, N. (2025). Diffusion-adsorption of water vapor in chemically
  activated lunar soil. *Planet. Sci. J.*, 6(7), 164.
  doi:10.3847/PSJ/ade5b2.
  Both Schörghofer papers are CC BY 4.0. Their Eq. 4 and Eq. 9 are used in
  U24-c.
- For the Do choice, the two originals:
  - Joens, J. A. (2001). The dissociation energy of OH(X²Π₃/₂) and the
    enthalpy of formation of OH(X²Π₃/₂), ClOH, and BrOH from thermochemical
    cycles. *J. Phys. Chem. A*, 105(49), 11041–11044. doi:10.1021/jp011833u.
    This is the source of 35584(10) cm⁻¹.
  - Sun, G., Zheng, X.-F., Qin, Y., Song, Y., Zhang, J., Amero, J. M., &
    Vázquez, G. J. (2020). Two-photon dissociation dynamics of hydroxyl
    radical. *Chin. J. Chem. Phys.*, 33(2), 129–134.
    doi:10.1063/1674-0068/cjcp2003026. This is the source of
    35580(15) cm⁻¹.
  - ATcT 35580 ± 3 cm⁻¹ as cited by Sun et al. 2020: Ruscic, B. (2015).
    Active Thermochemical Tables: sequential bond dissociation enthalpies
    of methane, ethane, and methanol and the related thermochemistry.
    *J. Phys. Chem. A*, 119(28), 7810–7837. doi:10.1021/acs.jpca.5b01346.
    Current values are at https://atct.anl.gov, which returned 403 to
    automated access.

References for the first batch:
- Tan, K. H., Brion, C. E., van der Leeuw, Ph. E., & van der Wiel, M. J. (1978).
  Absolute oscillator strengths (10–60 eV) for the photoabsorption,
  photoionisation and fragmentation of H₂O. *Chemical Physics*, 29(3),
  299–309. doi:10.1016/0301-0104(78)85080-0. The Crossref title field reads
  "H20s", which looks like a metadata typo.
- Kronebusch, P. L., & Berkowitz, J. (1976). Photodissociative ionization in
  the 21–41 eV region: O₂, N₂, CO, NO, CO₂, H₂O, NH₃ and CH₄. *International
  Journal of Mass Spectrometry and Ion Physics*, 22(3–4), 283–306.
  doi:10.1016/0020-7381(76)80088-5. This is the D04 branch source for
  H₂O04–06.
- Eland, J. H. D. (1975). The predissociations of water ions. *Chemical
  Physics*, 11(1), 41–47. doi:10.1016/0301-0104(75)80036-x.
- Powis, I., & Reynolds, D. J. (1991). Dissociation dynamics of D₂O⁺ ions
  prepared in the B̃ ²B₂ state. *Journal of the Chemical Society, Faraday
  Transactions*, 87(7), 921. doi:10.1039/FT9918700921. It studies **D₂O⁺**,
  so its use for H₂O falls under decision E4. Crossref gives only the first
  page.
- Banna, M. S., McQuaide, B. H., Malutzki, R., & Schmidt, V. (1986). The
  photoelectron spectrum of water in the 30 to 140 eV photon energy range.
  *Journal of Chemical Physics*, 84(9), 4739–4744. doi:10.1063/1.450008.
- Truesdale, C. M., Southworth, S., Kobrin, P. H., Lindle, D. W., Thornton,
  G., & Shirley, D. A. (1982). Photoelectron angular distributions of H₂O.
  *Journal of Chemical Physics*, 76(2), 860–865. doi:10.1063/1.443058.
- Liu, J., Salumbides, E. J., Hollenstein, U., Koelemeij, J. C. J., Eikema,
  K. S. E., Ubachs, W., & Merkt, F. (2009). Determination of the ionization
  and dissociation energies of the hydrogen molecule. *Journal of Chemical
  Physics*, 130(17), 174306. doi:10.1063/1.3120443. This is the para-H₂
  D₀ source needed for H₂O02 and H₂O05 (decision B2).

## 7. What this proposal does not claim

- It makes no change to D04, CHEMICAL_DATABASE, the roadmaps, test IDs or
  thresholds.
- No molecular or surface run was performed.
- Every proposed check is supporting evidence at most. Formal I22 also needs
  the complete family with approved kinetics, plus U15, U24 and
  same-campaign I13–I20.
- Observational comparisons (D07, D09, D15, LADEE water) are out of scope.
