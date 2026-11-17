# M3 H2O/OH photochemistry: evidence and implementation decisions

Research proposal, 2026-11-17. **Not approved physics, not an implementation,
and not a declaration that U15/U24/I22/I31 passed.** This work changes no
production source, observation, reference, acceptance limit or raw data.

## Recommended first scope

Proceed with the separately approved atomic M3 work. For molecules, first
implement the complete frozen **rate/channel inventory and nucleus/charge
bookkeeping**, with molecular transport disabled by default and rejected when
its required kinetic evidence is missing. A closed-box parent-survival/channel-
count check can verify these rates without inventing daughter velocities.
It cannot qualify daughter transport, escape, line brightness or water
migration. Full molecular I22/I31 and the combined U24 gate remain SKIPPED.

Do not enable the legacy mean-energy/isotropic daughter recipe as physically
qualified solar H2O/OH transport. Published work establishes wavelength- and
state-dependent product distributions; the frozen D04 files do not contain
those joint distributions. Replacing a distribution by its mean preserves
neither its speed distribution nor a nonlinear observable such as escape
probability. This is a missing-data/model contract, not a request to run longer.

## Authority, files and current scope

Checkout revision at inspection:
`17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`, with uncommitted development.
Roadmaps at the initial read differ:

| Document | SHA-256 |
|---|---|
| `srcMoon/development_and_validation.md` | `eba1571968afff4f948505b77578ee782c4fb54d02ac786ab9ecaa7ea00438d2` |
| `/data/vtenishe/moon_validation_data/development_and_validation.md` | `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4` |

The source copy contains the approved development-order/budget and later
implementation-status amendments absent from the external copy. Both retain
R3.2's species/rate/branch requirements and R4.7's later water migration scope.
The user's new M3 development/local-verification authorization controls work
order; it does not close M1/M2 or open formal observational qualification.
Native D04 labels and the external data AGENTS control data interpretation.

Read-only evidence is under
`/data/vtenishe/moon_validation_data/photochemistry/`. The package README,
provenance and QA were read. QA reports extraction/units/traceability PASS;
that is not qualification of a daughter kinetic distribution. The local
experimental OH product is `OH_exp`, native PHIDRATES identifier **`OHe`**.
Its physical parent is OH; it is neither a new isotope nor an approved
rovibrational-state population. Preserve that experimental provenance rather
than substituting another OH compilation.

Important frozen products:

| File in `processed/phidrates/` | SHA-256 |
|---|---|
| `photo_rates_1au.csv` | `3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847` |
| `reaction_channels.csv` | `2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061` |
| `cross_sections.csv` | `e900eb5514b616dcbbde7e3ca124b5c06d4d22a28a3cd429ed58a14c1933c860` |
| `branching_ratios.csv` | `430bb68558a43973f0622094049c87004dd522a6288a48935027856f48229ba7` |
| `thresholds.csv` | `a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1` |

Six processed PHIDRATES products and ten H2O/OH raw input descriptors were
independently hashed against the archived provenance/processing summary:
**16 matched, zero mismatched**. Numeric inventory uses the committed package
processor's products; no numerical value was digitized, OCRed, altered,
rejected or filled. Raw BrnOut headers and quiet-Sun text were also inspected.

The reproducible read-only inspection and its complete identities/rows are in
[`audit.json`](../../test_output/srcMoon/m3-molecular-research-20261117T072725.641599Z/audit.json),
beside `audit.py`. Reproduce this inventory from the repository root with:

```bash
/usr/bin/python3 -B test_output/srcMoon/m3-molecular-research-20261117T072725.641599Z/audit.py "$PWD" "$PWD/test_output/srcMoon/m3-molecular-research-20261117T072725.641599Z"
```

This regenerates an inspection artifact, not processed observations or a
reference solution. The only evaluated success was 16 archive identity matches;
no build, MPI launch or U/I test was run. The audit separately records file
hashes at its execution time because other approved development is concurrent.

## What the frozen rates actually specify

All following values are at 1 AU, optical depth zero. Quiet/active solar
endpoints are separate archived cases. `E` is the source's scalar mean excess
energy in eV, **not a sampled energy or angular distribution**.

| Parent and source products | Quiet rate, s^-1 | Active rate, s^-1 | Quiet/active E, eV |
|---|---:|---:|---:|
| H2O → H + OH | 1.03e-5 | 1.76e-5 | 3.41 / 4.04 |
| H2O → H2 + O(1D) | 5.97e-7 | 1.48e-6 | 3.84 / 3.94 |
| H2O → O + H + H | 7.54e-7 | 1.91e-6 | 0.699 / 0.697 |
| H2O → OH+ + H + e | 5.54e-8 | 1.51e-7 | 18.6 / 23.2 |
| H2O → O+ + H2 + e | 5.85e-9 | 2.21e-8 | 36.5 / 39.8 |
| H2O → H+ + OH + e | 1.31e-8 | 4.07e-8 | 25.0 / 30.5 |
| H2O → H2O+ + e | 3.31e-7 | 8.28e-7 | 12.4 / 15.2 |
| OH_exp → O(3P) + H | 1.20e-5 | 1.37e-5 | 2.00 / 2.14 |
| OH_exp → O(1S) + H | 8.33e-7 | 2.11e-6 | 10.0 / 10.0 |
| OH_exp → OH+ + e | 2.43e-7 | 6.43e-7 | 19.4 / 23.6 |
| OH_exp → O(1D) + H | 7.01e-6 | 1.76e-5 | 7.73 / 7.74 |

These are all seven/four archived reaction channels; none has been discarded.
The H2O three-neutral channel's source rate label is simply `O`, so its
electronic state requires the original channel reference, rather than an
assumed substitution. Several channel-specific reference blocks are empty,
even where the total-absorption block cites experiments. Keep those gaps.

For a network consisting of these **complete listed reactive families**, the
arithmetic channel-rate sums are:

| Family | Quiet sum, s^-1 | Active sum, s^-1 |
|---|---:|---:|
| H2O | 1.205635e-5 | 2.20318e-5 |
| OH_exp | 2.0086e-5 | 3.4053e-5 |

The CSV's `total_photoabsorption` cross section is not an additional reaction.
Do not add it to the partial hazards or assume every absorption destroys a
parent. Sampling `p_i=k_i/sum(k_i)` is the conditional selection of the complete
declared reaction family; it does not change any input rate. If some products
are not transported, retain their original hazards and explicitly account for
their products/reservoir transfers, or reject that requested model. Never
delete those hazards and rescale the easier transported channels.

The wavelength tables are rounded source values, not exact partitions at
floating-point precision. Among the 172 H2O and 64 OH_exp wavelength records,
the partial sum differs from `Total` in 104 and 47 records respectively. The
largest fractional differences found are -0.00294117647 at 114 and 121.5 nm,
and -0.002 at 120 and 121.5 nm. These arithmetic findings do not authorize renormalization or an
invented physical uncertainty. Preserve the numbers and distinguish source
printing precision from numerical summation error.

There is a further qualification issue: OH_exp's raw table has tiny nonzero
`O1S/H` values extending to 61320 Å (6132 nm, photon energy about 0.202 eV).
For a ground-state OH parent this cannot supply even the O(1S) excitation
energy, before paying dissociation energy. Such source tails/floors need an
authoritative interpretation; their physical role is not documented locally.
No clipping, threshold repair or silent extrapolation is proposed. The local
`thresholds.csv` supplies numerical entries for the four atoms, not the H2O/OH
channel thresholds required for an independent wavelength-resolved check.

## Energy and direction evidence

[Huebner & Mukherjee (2015), doi:10.1016/j.pss.2014.11.022](https://www.sciencedirect.com/science/article/abs/pii/S003206331400381X)
defines a spectrum/cross-section-weighted mean excess energy above a specified
channel threshold. The paper identifies photoelectron kinetic energy as the
principal ionization-energy recipient. A scalar average does not specify its
spectrum, a three-product energy split, or a joint distribution of molecular
internal excitation and translation. A threshold for a specified product state
must not have that same excitation energy subtracted a second time.

[Wu & Chen (1993), doi:10.1029/92JE03016](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/92JE03016)
construct solar H/OH velocity distributions using partial cross sections,
OH internal-state distributions and both solar continuum and bright lines.
Their published calculation reports multiple velocity groups and bounds where
inputs were missing. It is a useful independent benchmark candidate, not
permission to adopt one fixed speed or transfer its uncertainty to another
parent population/spectrum.

[Harich et al. (2000), doi:10.1063/1.1322059](https://doi.org/10.1063/1.1322059)
is a candidate state-resolved Lyman-alpha benchmark identified in the primary
Chang article's references. The publisher's full article was inaccessible
during this review; its numerical state/angle distributions have not been
acquired or frozen here. A monochromatic benchmark would retain its actual
parent-state and polarization conditions; it would not represent the full
solar endpoints. Its inaccessible original must not be treated as acquired
kinetic evidence.

[Chang et al. (2019), doi:10.1038/s41467-019-09176-z](https://pmc.ncbi.nlm.nih.gov/articles/PMC6424997/)
measures wavelength-dependent translational/internal-state structure, including
rotationally excited OH products. The primary article distinguishes parallel
and perpendicular laser-polarization measurements and gives a data-availability
route. This demonstrates why one isotropic, monoenergetic solar daughter law
cannot be inferred from a PHIDRATES mean. Its plotted spectra were not
digitized, and no distribution from that paper was implemented.

[van Dishoeck & Dalgarno (1984), *Photodissociation of OH*](https://adsabs.harvard.edu/pdf/1984ApJ...277..576V)
distinguishes molecular electronic pathways yielding ground or excited atomic
oxygen. [Nee & Lee (1984), doi:10.1063/1.447387](https://doi.org/10.1063/1.447387)
is the experimental absorption reference in the frozen OH total block.
An absorption measurement is not itself the missing daughter velocity/state
distribution. Neither source licenses filling the archived blank per-channel
reference fields.

For a future data-resolved implementation, the required object is, schematically,
`P(product internal states, translational energies, relative directions |`
`parent internal state, photon wavelength/polarization, channel)`.
The solar spectrum and partial cross section determine the photon/channel
sampling; the conditional distribution determines the daughter kinematics.
Both pieces require identified versions, units and independent oracles.

An isotropic kernel can remain a controlled mathematical verification case.
It is not automatically the physical angular law. For example, transforming
a measured electric-vector anisotropy to an unpolarized solar-beam angular law
requires an explicit polarization average; merely choosing isotropic directions
changes that law. No anisotropy coefficient or energy PDF is assigned here.

## Excited oxygen and OH relaxation

Do not collapse O(3P), O(1D) and O(1S) merely because their nuclear masses are
the same. They differ in internal energy, radiation and potentially chemistry.
A proposed representation is an explicit internal-state tag with separate
populations and decay events; distinct species IDs are another implementation
choice if their weights/movers and state conversion are qualified.

The primary official
[Ralchenko/NIST atomic-data presentation, ICTP/IAEA 2024](https://amdis.iaea.org/media/presentations/ICTP2024-Ralchenko-atomic-data-for-plasmas.pdf)
lists O(1S)→O(1D) at 5577 Å, smaller direct ground-term channels near 2972/2958 Å,
and O(1D) red/weak ground-term transitions. Summing its listed probabilities
illustrates lifetimes of roughly 0.75 s and 134 s. These are **not adopted
constants**: freeze an explicit complete
[NIST ASD](https://physics.nist.gov/asd) release/query and underlying references
before implementation, sum every independent radiative contribution once,
then use `tau=1/sum(A_ul)` and branch probabilities `A_ul/sum(A_ul)`.
O(1S) must not be forced to cascade through O(1D) with probability one.

Model radiative decay as a competing timed event, preserve the atom and its
trajectory, update its state, and account for photon energy/momentum. An
instantaneous-ground-state approximation needs a demonstrated space/time and
kinetics error budget; it is not implied by a short lifetime. Collisional
quenching is a separate process. Do not assume it negligible without the
selected density/rate bounds, or replace chemical quenching by spontaneous
decay. OH(A) radiative/predissociation and vibrational/rotational relaxation
likewise need an explicitly selected state model if those populations are
created; a ground-state OH solar table is not automatically valid for every
newly born OH state.

## Production code that motivates the decisions

The generic file
[`PhotolyticReactions_Huebner1992ASS.cpp`](../../src/models/photolytic_reactions/PhotolyticReactions_Huebner1992ASS.cpp)
contains the seven H2O rows and four OH rows, with compile-time quiet/active
selection. H2O's O(1D) row is represented as ordinary oxygen:

```cpp
{_OH_SPEC_,_H_SPEC_,-1},
{_O_SPEC_,_H2_SPEC_,-1},
{_O_SPEC_,_H_SPEC_,_H_SPEC_},
```

All three OH atomic-oxygen channels likewise share the same species ID:

```cpp
{_O_SPEC_,_H_SPEC_},
{_O_SPEC_,_H_SPEC_},
{_O_SPEC_,_H_SPEC_},
{_OH_PLUS_SPEC_,_ELECTRON_SPEC_}
```

The generic selector sets a channel's one scalar energy and calls a two- or
three-body generator. Its two-body direction is explicitly isotropic:

```cpp
if (ExcessEnergyTable!=NULL) TotalEnergy=ExcessEnergyTable[ReactionChannel];
// ...
double cosTheta=2.0*rnd()-1.0;
```

Those generators can verify mathematical energy/momentum closure for that
specified model. Their conservation does not identify the physical energy-
partition probability density. They also filter daughters whose compiled
species ID is negative; that does not provide a missing-species nucleus,
charge, energy or momentum ledger.

[`Exosphere_Chemistry.cpp`](../../src/models/exosphere/Exosphere_Chemistry.cpp)
creates each product species separately using its yield, weight and time step:

```cpp
double anpart=ProductionYieldTable[spec][specProduct]*
    (1.0-exp(-ProductTimeStep/ParentLifeTime))*ParentParticleWeight/ProductParticleWeight;
// each product loop samples its own channel containing that product
PhotolyticReactions::GenerateReactionProducts(spec,ReactionChannel,
    nReactionProducts,ReactionProductsList,ReactionProductVelocity);
```

Parent survival is then sampled separately. This may implement a statistical
yield model; it does not establish correlated daughters from one physical
event or event-by-event conservation. A new molecular adapter needs an
explicit common weighted reaction event or a separately justified statistical
conservation contract. Do not demand exact paired-event identities from this
existing independently sampled product algorithm without changing the contract.

The generic lifetime cache is initialized once per species; its default
`TotalLifeTime` does not use position-dependent illumination. The lunar
[`Moon.cpp`](../Moon.cpp) currently special-cases Na and delegates other
species to that generic processor. Its comment/code refuses unqualified Na+
product weights; it is not an established H2O/OH transport path.
[`Moon.h`](../Moon.h)'s existing lifetime callback also explicitly enables Na
only. `input/species.input` lists OH/H2O definitions, but that alone does not
compile or qualify them in the normal Moon application.

Future gas photon rates must use fluence/visibility at the **particle's**
actual position and epoch, with the selected spectrum and inverse-square
scaling. Surface projected incidence, absorptivity or a surface temperature
is not a gas photodestruction factor. Reuse the production geometry service
with the correct volume observation operator, not a duplicate test shadow law.

## Required bookkeeping and bounded verification

For every event preserve reaction identity, physical particle weight,
stoichiometric multiplicity, electronic/internal states, products not
transported, photon input and all reservoir transfers. H2O reactions conserve
two H nuclei and one O nucleus; OH reactions conserve one of each. Every
ionizing branch creates one positive ion and one electron. The two-H channel
has multiplicity two, not one species occurrence. A missing electron can be
an explicitly declared external reservoir; it cannot silently disappear from
charge/energy accounting.

Given an independently sampled two-body translational energy, conservation
fixes equal/opposite centre-of-mass momenta with
`p=sqrt(2*mu*E_trans)`; it does not fix the direction PDF or molecular internal
state. Three-body conservation similarly does not determine a physical
energy-sharing distribution. Incoming photon recoil, emitted photons,
photoelectrons, ion masses and internal/binding-energy conventions need a
consistent nonrelativistic bookkeeping contract. Any neglected recoil or
untracked carrier requires an explicit approximation/error budget.

Extend existing U15/U24/I22 meanings; assign no new permanent ID. A proposed
bounded first campaign has frozen endpoint/hash/topology checks, analytic
closed-box survival/branch counts, independent nucleus/charge identities,
missing-product refusal and weighted production-wiring diagnostics. Preserve
the roadmap's existing <=2% linked rate and 95% branching gates; U15's local
contract separately declares 99% branch-sampling intervals. These distinct
declared gates must both remain intact. Kinematic
tests later require measured/state-resolved distribution comparisons in
addition to energy/momentum invariants; invariants alone do not validate a PDF.
All pure fixtures must be labelled, separated from observations, and registered
through the existing individual/global runners. No actual test was executed
or implementation made by this research work.

## Decisions needed before molecular transport

1. **Approve the bounded staged scope above**, retaining full molecular
   transport as disabled/NOT QUALIFIED until kinetic inputs exist; or request
   a specifically identified laboratory/spectral benchmark first. A single-
   wavelength benchmark must not be renamed full solar transport.
2. Select the parent electronic/rovibrational population and authoritative
   channel-resolved energy/angular/internal-state files, including the three-
   neutral and dissociative-ionization branches. The cited papers are concrete
   acquisition candidates, not a complete frozen distribution dataset.
3. Select explicit excited-O/OH state tracking and a frozen radiative dataset,
   including all decay branches; define supported quenching/secondary chemistry
   or an independently bounded omission. No lifetime/model is approved here.
4. Choose a complete transported product set and explicit electron/photon/
   untransported-product reservoirs, with weighted event and recoil conventions.
5. Resolve the PHIDRATES molecular threshold/tail interpretation and blank
   per-channel provenance. Keep experimental OH identified as such and do not
   repair or renormalize the frozen evidence to make a test pass.

No approval is inferred by this document. Missing answers/data must keep the
affected capability/test SKIPPED or disabled while independently authorized
atomic M3 implementation and local verification proceed.
