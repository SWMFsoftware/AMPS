# M3 ion-surface and charge-exchange kinetics: data-source proposal

Research date: 2026-11-17. **PROPOSED acquisition/qualification choices only.**
No reaction execution, surface law, observations, source-data readiness, or
acceptance criteria were changed. Downloaded numerical arrays and papers in
the acquisition packet remain research evidence; they are not a READY Dxx package.
Canonical `/data/vtenishe/moon_validation_data/**` was read-only throughout.

The companion [atomic and molecular kinematics source proposal](M3_KINEMATIC_REFERENCES_PROPOSAL.md)
covers Na photon/electron recoil and H2O/OH daughter distributions. This proposal
covers the separate ion-surface and charge-exchange data requirements.

The durable acquisition packet is
[`test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z/`](../../test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z/).
Its original downloads, retrievers and manifests are retained without modifying
canonical validation data. Source inspection revision is
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`, with concurrent approved M3 work.

## Recommendation for the next decision

Keep the implemented He/Na **total-rate** CEX registry separate from full
particle CEX. For the next surface research package, prioritize He+ on silica
and silicate analogues, with Apollo-sample microscopy as a separate reservoir
constraint. These are the strongest concrete sources found here; they do not
justify a universal He/Ne/Ar/Na ion sticking coefficient. Obtain the missing
joint state/angular CEX data and incident-ion outcome data before enabling a
full kinetic operator. No preferred source supplied by the user was assumed.

| Proposed package | First authoritative source | What it can support | What stays unresolved |
|---|---|---|---|
| He+ charge survival on a declared SiO2 target | Téllez et al.2012, Table I and measured LEIS spectra | A material/geometry/energy-specific charge-survival submodel, after exact numerical transcription and QA | Total reflection vs implantation, arbitrary emission angle, real lunar mixtures, and low-energy extrapolation |
| He+ scattering on a silicate target | Brötzner et al.2026, CaSiO3, 1–3keV | Charge fractions resolved by scattering class, after acquiring authors' numerical spectra and documented analysis | The authors state numerical supporting data are request-only; no request was sent |
| He implantation/microstructure on altered enstatite | Cymes et al.2026 plus Zenodo15978325 | Named-target morphology/composition/spectroscopy constraints and an independently labelled BCA comparison | Incident-ion trapping yield, outgoing charge/energy/angle probabilities, general residence law |
| Actual lunar He reservoir | Cymes et al.2024 plus Zenodo10794539 | Evidence for trapped helium in specified Apollo17 grains and native EELS analysis | Event-level implantation probability and time-dependent thermal-release coefficients |
| H+He CEX angular benchmark | Johnson et al.1989 plus original Rice table | A measured **5keV**, finite-angle, state-summed benchmark | All-angle closure and joint electronic-state/recoil distribution |
| H+He CEX electronic channels | OPEN-ADAS ADF01 He donor/H receiver data | Explicit n/l/m-resolved integral source data over its native 70–15000keV/u domain | Solar-wind-keV coverage, a differential angular law, and independently normalized joint probabilities |
| H+Na CEX joint kinetics lead | Courbin-Gaussorgues et al.1990 | State-to-state differential theory over 0.5–5keV, if actual numerical arrays become available | Full publication/arrays were not acquired; cannot infer a sampler from its abstract |

## 1. Ion-surface sources and their limits

### 1.1 First numerical charge-survival source: silica / He+

[Téllez et al.2012, APL101151602, DOI10.1063/1.4758699](https://doi.org/10.1063/1.4758699)
has an accessible [author-institution version of record](https://pure.tue.nl/ws/portalfiles/portal/3503988/583864115491353.pdf).
It measures He+ on isotope-exchanged silica at0.6–7keV, normal incidence and
145° scattering. Table I gives numerical characteristic velocities and fit
ranges for Si/O and distinct neutralization regimes. The intermediate regime
must not be filled by assuming either limiting law. The observed charge
fraction conditions on a detected scattered particle: it is not the fraction
of all incident ions that implants, sticks, or returns as a neutral.

**Acquired:** original5-page PDF,754543bytes,
SHA256`bf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2`.
Table values have been inspected in text, but no machine-readable production
coefficient file or independent numerical oracle was created in this task.
Proposed use: qualify a declared silica/He charge-survival submodel first,
with its actual angular acceptance and source regimes retained. It cannot
be applied to an arbitrary lunar regolith facet without a reviewed material
and roughness/geometry model.

### 1.2 More relevant silicate scattering; numerical data not public

[Brötzner et al.2026, PRB113075405, DOI10.1103/3ks5-1c83](https://journals.aps.org/prb/abstract/10.1103/3ks5-1c83)
uses1–3keV He+ on CaSiO3 with a145° detector. It separates single, dual and
multiple scattering using a dynamic binary-collision model, with measured
spectra constraining charge fractions. The source explicitly says supporting
numerical data are available from the authors on request. The extracted
fractions are conditional on that experiment and analysis, not a measured
universal implantation yield. Simulation output must remain distinguishable
from measurements.

**Acquired:** [institutional publisher PDF](https://repositum.tuwien.at/bitstream/20.500.12708/226480/1/Broetzner-2026-Physical%20Review%20B-vor.pdf),
3156406bytes,SHA256`11a94ab8d81f99329fa6423f3adf41e892bcce9f82fb98f684b406ac90518bb6`.
**Not acquired:** authors' numerical spectra/fit covariance/response/BCA input
decks. Requesting those files is the concrete next acquisition route; no
external communication was authorized or sent.

### 1.3 Enstatite implantation and space-weathering data: useful but different observables

[Cymes et al.2026, PSJ7,6, DOI10.3847/PSJ/ae2657](https://doi.org/10.3847/PSJ/ae2657)
has a public author-deposited dataset
[Zenodo15978325,v1](https://zenodo.org/records/15978325), also linked by
[NASA NTRS20250006770](https://ntrs.nasa.gov/citations/20250006770).
The deposit explicitly labels microscopy, EELS/EDS CSVs and SDTrimSP text.
Those are distinct data classes. Its depth/composition text is a simulation,
not independently measured ion trapping. Raw EELS spectra alone do not
determine absolute incident-ion retention without the experimental dose,
geometry, sensitivity, selected regions and processing protocol.

**Acquired:** exact Zenodo metadata and three small original files; all three
match archive MD5 and byte counts. `enstatite_selected_0.txt` is a native
depth[nm]/elemental profile; the two `selected_1/2.csv` files are raw He/matrix
spectra with no self-contained column labels. Their axes must be recovered
from authoritative article/instrument metadata before processing. The
publisher HTML request returned a bot challenge, not the article, which is
recorded rather than treated as acquired scientific text. Most of the467MB
deposit was deliberately not downloaded for this bounded review.

### 1.4 Lunar reservoir constraints and heating: no universal residence law

[Cymes et al.2024, DOI10.1038/s43247-024-01349-z](https://www.nature.com/articles/s43247-024-01349-z)
has a [NASA-hosted accepted manuscript](https://ntrs.nasa.gov/citations/20240003651)
and [original Apollo17 EELS deposit, DOI10.5281/zenodo.10794539](https://zenodo.org/records/10794539).
The experiment identifies He in vesicular metallic-iron particles in specified
lunar grains. It constrains reservoir/material interpretation; it does not
provide a complete incident-energy/angle dependent implantation law or a
unique thermal-release timescale. The accepted manuscript and archive
metadata were acquired; the682.5MB source archive was not downloaded.

[Kling et al.2025, NASA NTRS20250000214](https://ntrs.nasa.gov/citations/20250000214)
is an extended conference abstract on heating a specified Apollo ilmenite
sample. It is useful evidence of distinct trapping/release behavior, not a
ready coefficient table. The two-page original was acquired. Its cited
exospheric timescale must not be adopted as an all-material microscopic
residence time. [UWFDM1095](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1095.pdf)
was also acquired as a historical primary technical-report lead; its scanned
body was not numerically transcribed, so no diffusion coefficients are
claimed here.

No numerical Na+/Ne+/Ar+ on lunar-mixture **joint** outcome kernel was found
and qualified in this bounded search. That statement is an acquisition gap,
not a claim that the literature contains none. Sodium migration inside an
irradiated semiconductor oxide, neutral-Na adsorption, ion-induced sputtering
and incoming Na+ charge survival are different experiments and cannot be
substituted for one another.

## 2. CEX product-state/differential sources

### 2.1 Acquired primary angular data near keV energies

[Johnson et al.1989, PRA403626, DOI10.1103/PhysRevA.40.3626](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.40.3626)
reports proton-He charge-transfer scattering at5keV. The originating
[Rice group data index](https://www.ruf.rice.edu/~atmol/charge_data.html)
supplies the [numerical source table](https://www.ruf.rice.edu/~atmol/data/charge_transfer/H%2B/H%2B_He.txt).
Its accompanying index establishes laboratory angle and10^-16cm2sr^-1 units.
**Acquired:**1948bytes,SHA256`fe022c392c86c7bdb8e90c40c93286963ed77e3ef6466712e237dc7e74d577e3`.
The raw table contains distinct angle uncertainties and duplicated-angle
records; these must remain intact. Its first tabulated angle0.002° differs
from the article abstract's stated0.02° lower bound. Resolve that discrepancy
from the original paper/author metadata before qualifying an integration
range. Neither missing angular tails nor unresolved states may be silently
renormalized into a complete collision sampler.

The same primary index provides proton-Ar angle tables at0.5/1.5/5keV;
the7618-byte original was acquired. It is a concrete **future** Ar differential
source, not approval to add an Ar total-rate law or assume electronic-state
closure. Native blanks are preserved.

### 2.2 Acquired state-resolved integral data; wrong range for ordinary keV wind

[OPEN-ADAS ADF01 specification](https://open.adas.ac.uk/adf01) describes
state-resolved **integral**, rather than scattering-angle differential, cross
sections. The actual
[He(ground) donor/H+ receiver record](https://open.adas.ac.uk/element/detail/adf01/qcx%5D%5Bhe0/qcx%5D%5Bhe0_gyt%5D%5Bh1.dat)
and its original file were acquired. The file contains27 tabulated energies
from70 to15000keV/u, n1–5 totals, nl/nlm partials and explicitly documented
high-n/asymptotic completion. Those hierarchical sums are not independent
branches to add together. Numerical state coverage does not supply angular
correlations. No high-energy data or asymptotic recipe should be extrapolated
to a~keV solar-wind calculation.

Original file:18328bytes,
SHA256`b97eb243b62ddf124cab9782b1247090ce819afe08abd4050d3262b195e3ad4c`.
**Proposed:** retain as a high-energy reference/lead, seek a source-supported
low-energy joint dataset, and keep particle CEX OFF in unsupported regimes.

### 2.3 State-selective differential literature leads

[Schöffler et al.2007, arXiv0711.4920](https://arxiv.org/abs/0711.4920)
provides experimental state-selective differential capture for proton/helium
projectiles at60–630keV/u. The original17-page PDF was acquired
(SHA256`f92aa62997fb692018f8c0b935245bde84cd9a149a8546a1a985bfc50c99081e`).
It is a useful real joint-kinetics benchmark at those energies, not an
ordinary-keV replacement. Complete machine-readable joint arrays were not
acquired.

For Na, [Courbin-Gaussorgues et al.1990, DOI10.1088/0953-4075/23/21/023](https://doi.org/10.1088/0953-4075/23/21/023)
and the [author-institution record](https://oa.upm.es/61369/)
identify state-to-state differential **calculations** for ground/excited Na
over0.5–5keV. The latter download returned a bot challenge, not the paper;
the primary abstract remains a lead. Obtain actual numerical arrays and state
definitions before assessing joint closure. Ground/excited-target mixtures
must be explicit; a3p-selected target is not a3s default.

[Jaecks et al.1970, DOI10.1103/PhysRevLett.25.491](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.25.491)
measures the H(2p) capture component at4–20keV in proton-He collisions. This
partial-state benchmark does not complete all channels. It was identified
from the primary publisher, but full numerical arrays were not acquired.

## 3. Required contracts before implementation activation

**Ion-surface event:** distinguish incoming charge neutralization, escaping
ion/neutral reflection, implantation/trapping depth, subsequent diffusion or
release, and sputtered target material. A charge-survival probability
conditions on a scattering history; a sticking coefficient cannot represent
all these processes. A production contract needs incoming species/charge,
energy/angle domain, target composition/structure/temperature/irradiation,
joint outgoing outcomes, and surface charge/electron/energy ledgers. No
isotope mixture, all-neutralization rule, albedo or residence constant is
inferred from the available papers. An idealized perfectly absorbing
collector would require an explicitly approved bounded experimental model;
it would not count as validation of lunar-regolith ion trapping.

**CEX event:** preserve projectile and donor identities; use state-dependent
energy defects and joint recoil/angle/internal-state distributions. A total
cross section, independent state partials and a separate angular marginal do
not uniquely specify their correlations. No factorization, forward-only
velocity swap, all-ground-state selection or unmeasured angular-tail closure
is approved here. A Maxwellian bath additionally requires an approved full
integration and treatment of its out-of-source-range tails.

The new production `ChargeExchange::Registry::GetMetadata` explicitly returns:

```cpp
bool daughterKineticsAvailable = false;
bool productionParticleExecutionEnabled = false;
```

That accurately describes the current source-supported scalar library. Its
`ColdBeamRateS(id, energyEvPerU, relativeSpeedMps, densityM3, ...)` accepts
explicit source-domain inputs and does not choose missing product physics.

## 4. Concrete proposed next work and decisions

1. **Recommend approving acquisition/QA first:** preserve the small source
   tables already acquired, document native variables/uncertainties and
   source inconsistencies, transcribe only published numerical tables with
   committed reproducible scripts, and obtain request-only spectra/joint
   arrays where necessary. Requesting authors' data requires separate explicit
   authorization; none was sent during this work.
2. **Recommend a He-first laboratory scope:** qualify a named He+/silica or
   CaSiO3 scattering/charge-survival operator within the actual source geometry
   and energy range. Keep implantation and lunar material extrapolation
   separate until supported by additional data. Na+/Ne+/Ar+ should remain OFF
   rather than inherit He parameters.
3. **Recommend retaining CEX rate-only scope for now:** integrate no full
   daughter executor until a reviewed all-channel/state/angle model is
   available. A separately approved reduced neutral-loss model could use the
   current total-rate laws, but must declare that it does not predict ENAs,
   pickup-ion velocities, or transported product populations.

These choices can be decided independently. None follows automatically from
approval of scalar rates or of the existing interpolation law.

## Evidence and qualification status

`retrieve.py`, `retrieve_more.py`, `retrieval.json` and
`retrieval_more.json` record exact URLs, final URLs, bytes, hashes and bounded
downloads. Sandbox DNS failures were retained separately before the authorized
host-network retrieval. Bot-response HTML is retained as acquisition failure
evidence, never scientific content. The three small Zenodo science files
match their producer MD5/size; that source-identity check is not validation.

Roadmap identities remain different:

| Document | SHA256 |
|---|---|
| Source roadmap | `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` |
| Canonical-data roadmap | `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4` |

Source amendments control approved early M3 development. Native source
metadata and common scientific/predecessor gates still control qualification;
the differences were not reconciled silently. Full lunar ion-surface and CEX
product validation remain **SKIPPED / NOT VALIDATED**. No application rebuild,
MPI simulation, observational scoring or automatic data promotion occurred
in this research task.


## Selected numerical-file identity audit

Producer MD5 verifies the selected native Zenodo files; SHA256 additionally binds
AMPS research evidence. All source names and URLs are in `retrieval_more.json`.
No science value was edited, regridded, interpolated or renormalized.

| Native file | Producer MD5 (verified) | Acquired SHA256 |
|---|---|---|
| `2b_He_Enstatite_Elements_Depth_Profile_SDTrimSP_compound.txt` | `47dbd21f2481968ed341008a36149cc2` | `6e2dc62c5d7b7f977d08184952b2c28b4af047c7fe692d401f49405967395078` |
| `9c_raw helium spectrum ROI9.csv` | `00514102a8f2653e350d25ac8b5c69c6` | `85ea998b5061861f2f78060047491f2ec3c7edbc7e78e4b3da6e6e99a112e720` |
| `9c_raw helium spectrum ROI9_matrix.csv` | `da0873c5205ab1c8a77eb7286f679e48` | `5bc6708e02c03fb7ef59e7ccebd0eff11e11ea9fa3dcd2d9ff9b8e4da7d75c49` |

The exact source-download scripts are in the acquisition packet. They use only
primary URLs, cap each response at20MiB, preserve failure responses, and never
write under the canonical data root. Re-running them is a research acquisition,
not automatic review or promotion. Recheck producer versions and science
metadata before generating any future processed package.
