# M3 physical evidence and community-model options — 2026-11-19

This follow-up prioritizes laboratory measurements and spacecraft observations,
then independently calculated atomic/molecular data and explicitly published
reduced models. It follows the user's request to find a practical, physically
documented route through the remaining M3 gaps. Source revision:
`8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129`, with current working changes.
No production law, frozen source selection, observational membership, acceptance
criterion, physical profile or installed chemical release is changed here.

**Some earlier “missing data” entries now have concrete solutions.** A physical
Na+ mass can be derived from measured source quantities. An actual community
O I transition matrix is available. Original water-fragment spectra and Apollo
water-desorption arrays have been acquired. Other gaps can be addressed by a
clearly named reduced model rather than requiring an unavailable complete
microscopic distribution for every observable. That model still needs its own
source, conservation, domain and approximation-error contract.

## 1. State, mass and atomic photon data

### Na/Na+ mass and the separate threshold crosswalk

The physical ground-state ion mass follows the measured neutral isotope mass,
electron mass and ionization energy:

```text
m(23Na+) = m(23Na) - m_e + I_Na/c^2
```

[CODATA 2022](https://physics.nist.gov/cuu/Constants/)
provides the consistent constants/convention;
[NIST isotope masses](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Na)
and the [Na ionization-energy source](https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=Na&units=1&format=1)
provide the measured inputs. The retained arithmetic gives nominal
`m(23Na+) = 3.8174499149007744e-26 kg`. Its first-order propagated standard
uncertainty must include covariance; the report supplies a conservative
triangle bound when covariance is unknown, rather than assuming independence.
This is a source-derived quantity, not an invented ion mass or a new measurement.
The original equal ExoAtom neutral/ion mass field is not used as physical truth.

The remaining Na issue is a **source convention and implementation decision**.
D04's authenticated threshold, `2412.58 Å`, has the retained processed value
`5.13907097104346 eV`; its difference from the selected modern ionization
energy is `5.98895654e-6 eV`. The rate HTML contains rate and mean excess energy,
not an independent `5.13907 eV` threshold. That truncated value in the earlier
atomic-source note was a display rendering; the source-tree note is corrected
here. The immutable r2 documentation retains its historical version. Propose an
explicit versioned crosswalk: retain the approved D04 hazard, preserve its
source-effective photon energy, use the selected physical ionization energy
in the event Q, and account for the corresponding change in available kinetic
energy. Do not reinterpret printed digits as an uncertainty, edit D04, or widen
the numerical guard. The effect on a newly recomputed photo-rate requires a
spectral/cross-section convolution; the small energy offset alone does not bound
that rate effect. Source permission/profile admission is a separate operation.

For photon angles, [Hosea et al.2023](https://doi.org/10.3390/atoms11100125)
discusses Na dipole/non-dipole behavior and the Cooper minimum. The existing
approved isotropic source-mean closure can remain a named reduced experiment.
A dipole `beta=2` Na3s case is a useful separate sensitivity model within its
domain, not a universal measured law over all wavelengths. Any extension must
retain electron/photon momentum and energy ledgers and independent model-error
bounds for the advertised observable.

### Oxygen cascade: an acquired community source with observational checks

[CHIANTI11.0.1 O I](https://db.chiantidatabase.org/o/o_1.html)
provides native `o_1.elvlc` and `o_1.wgfa`. The A values identify
Froese Fischer & Tachiev2004; the energy source identifies NIST2010. The acquired
five-lowest-level model contains the three O(1D) exits to the ground fine
structure and the three O(1S) exits, including weak channels often omitted.
Use `tau_u=1/sum_l A_ul` and `p_ul=A_ul/sum_l A_ul` for that declared model.
Its wavelengths are vacuum Å; do not relabel them as air wavelengths.

The selected matrix predicts an approximately16.07 photon-number ratio for the
green/297.2nm branches. Independent Mars observations report
[16.5±0.4, Gérard et al.2020](https://doi.org/10.1038/s41550-020-1123-2)
and [15.8±2.4, Soret et al.2023](https://doi.org/10.1029/2023JE007762), with different
uncertainty definitions. These favor a modern complete source choice over an
unexamined historical ratio near9.8. Do not mix a measured ratio with an unrelated
lifetime, equate photon and energy ratios, or invent independent errors.
Quenching/alignment and additional omitted-process bounds remain conditional
on the selected environment and observable.

This is not a claim that the matrix passes every independent line ratio.
[Sharpee & Slanger2006](https://doi.org/10.1021/jp056163x) reports a red-doublet
ratio of2.997±0.016, while the native CHIANTI A6300/A6364 ratio is approximately
3.09519. Its photon-versus-energy convention and the coherent alternative
transition calculations must be checked explicitly. Do not splice individual
A values or select a branch after examining AMPS residuals.

## 2. H2O/OH photochemistry: measured spectra and published approximations

The [Mordaunt et al.1994 Lyman-alpha experiment](https://doi.org/10.1063/1.466880),
also evaluated in [NASA/JPL19-5](https://jpldataeval.jpl.nasa.gov/pdf/NASA-JPL%20Evaluation%2019-5.pdf),
distinguishes OH(X)+H, O(1D)+H2, OH(A)+H and O(3P)+2H products. This supplies a
monochromatic state-resolved benchmark; it does not specify broadband solar
branching for every parent population.

Two useful original numerical archives were actually downloaded:

* [Chang et al.2021, three-body dissociation](https://doi.org/10.1038/s41467-021-22824-7):
  publisher Source Data contains `f1.xlsx`, `f2.xlsx`, `f4.xlsx`. It separates
  measured polarization/energy spectra, fitted component spectra and mixed
  input spectral curves. Nine92–109nm settings constrain O(3P)/O(1D)+2H and
  binary products; sequential behavior challenges a universal flat three-body
  energy split. Some source spectral units and signed wavelength entries need
  interpretation before integration. No fitted component is relabelled a raw
  measurement or converted into a complete joint probability distribution.
* [Chang et al.2021, vibrational H2 products](https://doi.org/10.1038/s41467-021-26599-9):
  publisher Source Data contains `Fig3.xlsx`/`Fig4.xlsx`, giving kinetic-energy
  spectra and fitted vibrational populations at six100–112.81nm settings.
  Substantial energy enters H2 internal excitation. This O(1S)+H2 channel is
  outside the frozen seven-channel D04 H2O family; adding it needs explicit
  source/network/test-contract review, not a silent extra hazard.

[Tabone et al.2021](https://doi.org/10.1051/0004-6361/202039549)
offers a concrete quantum-state/cascade model tested against Spitzer OH emission.
Its short-wavelength exclusions and sublevel assumptions are explicit. Newer
[JWST observations and author arrays, Zannese et al.2024](https://zenodo.org/records/10634184)
distinguish alternative OH sublevel-population models. The acquired Fig2 file
has2851 six-column rows: columns2–4 are observed/processed spectra; columns5–6
are different synthetic models. Units are native micron and MJy/sr; there is no
uncertainty column. Do not use its model columns as observations or compare a
lunar simulation directly to a protostellar environment.

There is also a documented coarse option: the comet
[vectorial/MonteCarlo literature](https://arxiv.org/abs/astro-ph/0702193)
uses an approximately1km/s OH recoil; [Wu & Chen1993](https://doi.org/10.1029/92JE03016)
instead derives multi-component solar H/OH distributions with internal-energy
accounting and bounded missing inputs. Fixed isotropic recoil is therefore a
published reduced model, not a newly measured universal distribution. For
lunar escape and state-sensitive brightness, compare the chosen approximation
against the actual fragment spectra before claiming adequate accuracy.

For resolved two-body events, subtract the selected internal-state/formation
energy once, then share translational energy by mass with opposite center-of-
mass momenta. For a randomly oriented parent in the electric-dipole limit, a
linear-polarization `1+beta*P2(cos theta)` law must be averaged over transverse
polarizations before representing an unpolarized solar beam. That analytical
operation does not generally give isotropy. Its assumptions and axis convention
must be declared; no beta value or missing wavelength support is invented here.

[Harper et al.2019](https://doi.org/10.1063/1.5091966)
supplies a separate OH photoionization measurement,12.6–15eV, with an absolute
13.8eV point and stated uncertainty. Its disagreement with earlier measurements
is real source uncertainty, not a reason to tune a cross section to AMPS.
[Suárez, Méndez & Rabadán2018](https://doi.org/10.1039/C8CP03725F)
provides state-dependent H2O+ fragmentation calculations checked against
coincidence experiments. Neither source licenses borrowing neutral dissociation
kinematics for every molecular-ion channel.

## 3. Electron ionization and charge exchange

The detailed [collision evidence report](../../test_output/srcMoon/m3-physics-literature-20261119T072531Z/collisions/report.md)
records primary measurements, native files and uncertainty conventions.

* [Opal et al.1971/72](https://doi.org/10.1063/1.1676707)
  gives measured secondary-electron distributions. A truncated Lorentzian-like
  energy-sharing law is used by published degradation/plasma models. Its width
  is species/source dependent; the literature uses different He values and no
  supported Na width was found. Isotropic/forward electron angles and a BED/BEB
  energy marginal remain separate approximations. Integral Voronov coefficients
  do not select those distributions, and a mixed rate/shape model needs review.
* [Murray2005](https://doi.org/10.1103/PhysRevA.72.062711)
  supplies measured Na electron-coincidence structure in restricted geometry.
  This can test a proposed differential model where measured; it is not a full
  solar-bath joint table.
* [Berkowitz & Zorn1984](https://doi.org/10.1103/PhysRevA.29.611)
  measures the H(2s) partial of proton–Na capture at500–2500eV. It adds useful
  lower-energy evidence but is neither total CEX nor ground-H production.
* [Crandall1970, NASA-CR-114270](https://ntrs.nasa.gov/citations/19700032941)
  contains actual differential/capture/state-partial numerical tables at3–20keV
  and1°–7° lab angles. The acquired PDF has real blanks and per-table errors;
  OCR/text aids are not an approved numerical transcription.
  [Johnson et al.1989](https://doi.org/10.1103/PhysRevA.40.3626)
  separately measures He CEX at5keV; its lower-energy CEX curves are theoretical.

The familiar nuclear-velocity-preserving ENA approximation is useful for
appropriate resonant/high-energy channels; see
[Lindsay & Stebbings2005](https://doi.org/10.1029/2005JA011298).
Na/He channels are not symmetric H+/H resonance. A forward two-heavy-product
model must still conserve total momentum and use `K_final=K_initial+Q` in the
center-of-mass frame. Copying both velocities unchanged when `Q!=0` generally
fails energy accounting. Collision partners must be sampled with the actual
`f(u)*sigma(g)*g` weighting or eliminated through an explicitly defined reduced
reservoir model; no finite table silently defines an unbounded bath tail.

### A practical alternative for heavy-ion observables

[Odelstad et al.2018](https://doi.org/10.1029/2018JA025542)
describes ions initially moving with their parent neutrals.
[Poppe, Halekas & Harada2022](https://doi.org/10.1029/2022JE007422)
uses an even coarser initial-velocity approximation for accelerated lunar pickup
ions. This supports proposing a separately named **inherit-parent-velocity,
unresolved-electron-reservoir** model for heavy-ion density/pickup observables.

Create the ion with the actual parent velocity/time/represented number and
resolved state. Transfer the exact signed heavy-particle charge, momentum and
energy differences to the electron reservoir, using the declared physical/reference
mass and formation-energy convention. Electrons remain unresolved and their
spectra/heating partition are not predicted. CEX requires a reservoir containing
the consumed plasma ion and its neutral/state products, not an electron-only bath.

This proposal does **not** require a full outgoing-electron joint PDF for its
limited observable. It does require a source-supported bound or an independent
higher-resolution comparison for omitted ion recoil and bath approximations.
Weak fields, low-energy bound trajectories, escape thresholds and light He can
be sensitive. Exact ledgers verify reduced bookkeeping, not microscopic recoil
accuracy. This variant is a new model/specification choice for review; the
existing strict source/provider guards are not removed by this research.

## 4. Surface response and water residence

[CENA measurements](https://doi.org/10.1029/2011JE004019)
support an empirical H energetic-neutral energy spectrum, with separate
[incidence-dependent angular measurements](https://doi.org/10.1029/2011GL049362).
The reported median neutral fraction is0.19, with quartiles0.16/0.21. These are
H and instrument/domain-specific results, not a Na/He/Ne/Ar reflection law or
a measurement of permanent retention by subtraction. Positive-ion reflection
is a different observable.

Binary-collision models such as SDTrimSP/TRIM are community tools for stopping,
reflection and implantation. [Szabo et al.2022](https://doi.org/10.1029/2022GL101232)
explicitly models porous enstatite; it does not calculate neutralization.
Measured sputter yields, outgoing charge fractions and incoming retained
fractions must retain separate meanings. An ideal collector remains appropriate
for mathematical first-impact verification, with no physical regolith claim.

The newly acquired [Jones et al.2020 author deposit](https://zenodo.org/records/3822365)
contains all14 native files,656794bytes, with producer checksums and CC-BY4.0
evidence. It includes actual Apollo14163/10084 TPD, fitted/inverted quantities
and simulations. This permits a real specimen/heating/coverage benchmark.
The inversion's prefactor and sticking-at180K assumption remain attached;
unit sticking and complete thermal accommodation are not measurements extracted
from these data. Native header ambiguity, signs and error bounds remain intact.

Two published reduced choices are available: a
[Jones et al.2018 water-cycle model](https://doi.org/10.1029/2018GL080008)
with unity sticking and a fixed55kJ/mol desorption energy/prefactor1e13s^-1;
and a [Schörghofer2023 coverage/site-energy model](https://doi.org/10.3847/PSJ/acf19b)
with an explicitly idealized efficient-diffusion end member. Prefer a named
specimen-derived site distribution for the laboratory benchmark; retain these
published alternatives for model sensitivity. Do not splice a prefactor from
one model into an energy inversion based on another. OH chemisorption/
recombination and molecular-H2O desorption remain different processes.

## 5. Measured plasma and independent observational products

The prior empty merged D06 package is not proof that electron data are absent.
Official [THB_L2_ESA metadata](https://spase-metadata.org/NASA/NumericalData/THEMIS/B/ESA/PT3S)
provides electron density `thb_peef_density` in cm^-3 and average temperature
`thb_peef_avgtemp` in eV, together with mode/cadence and validity information.
[THC_L2_MOM](https://spase-metadata.org/NASA/NumericalData/THEMIS/C/MOM/PT3S.html)
provides onboard electron moments as an alternative product. Actual CDF metadata,
quality flags, spacecraft-potential treatment and the distribution/model choice
must control processing. Measured moments are not proof of a Maxwellian shape.
Canonical campaign/probe selections remain empty; no arbitrary event is chosen
or imported by this research.

[Poppe2022's author deposit](https://zenodo.org/records/6994648)
contains model results and reduced ARTEMIS observations. The two small original
files were acquired: an observational non-solar-wind ion grid with distinct
all-ion/E-aligned/2–25keV flux columns, and a separate species-rate model table.
Fluxes are in cm^-2 s^-1, density in cm^-3 and spatial axes in source lunar radii.
Native zeros have no newly assigned missing/non-detection meaning. This is an
independent pickup-observable candidate after its native roles, selection/operator
and predecessor gates are qualified. Composition ambiguity prevents treating
an unresolved ion flux as a measured sodium-only distribution.

For future variable photon forcing,
[Moore & Mapaye2025](https://doi.org/10.3847/PSJ/add6a2)
provides TIMED/SEE-driven rates and published F10.7P fits for48reactions/12species.
Those are observation-driven rate calculations, not direct daughter distributions.
Their source inventory does not automatically cover Na/He/Ne/Ar or replace the
approved D04 quiet/active endpoints.

## 6. Recommended next work and exact decision boundary

1. Prepare an explicit Na source/state/constant crosswalk and compatible generated
   mass/charge configuration. Implement the energy-convention correction only
   after reviewing that source choice; retain canonical D04 unchanged.
2. Use the acquired CHIANTI oxygen matrix as a declared source model, preserving
   all listed exits and native units. Independently check sums, state mappings,
   decay budgets and modern observational photon ratios before claiming closure.
3. Propose the heavy-ion/reservoir variant for number/pickup observables and a
   separate source-resolved differential mode for recoil/electron observables.
   Freeze domains, source/model discrepancy and comparison observables before runs.
4. Prepare maintained importers for the newly acquired molecular/APollo numerical
   originals, preserving measured/fit/model roles and unknown fields. Publish a
   new explicit database release after review; never edit r1/r2 in place.
5. Benchmark named molecular/CEX/surface approximations against measurements in
   their actual wavelength, parent, material, energy and angle conditions. Keep
   full frozen hazard membership, source zeros/gaps and existing numerical gates.
6. Acquire appropriate electron CDFs once the scientific campaign/probe intervals
   are declared; qualify the real forcing/operator before event claims. Continue
   local prescribed-bath verification independently.

The existing statistical/numerical protocol and P01–P09 approvals stand. Newly
proposed source crosswalks and reduced-provider variants are not automatically
approved or implemented. Physical source/model uncertainties are distinct from
the1–3% implementation/convergence criteria; neither is widened to hide a failure.
Full I21–I24 still require actual linked evidence and original predecessors.

## 7. Evidence, reproducibility and status

Research outputs are under
`test_output/srcMoon/m3-physics-literature-20261119T072531Z/`:
`atomic/`, `collisions/`, `surface/`, `molecular/`, explicit parent acquisition
plans, actual original archives, source ledgers and inspection/completion JSON.
The maintained acquisition utility is
`utility/chemical_database/acquire_sources.py`; it records HTTP outcomes, UTC,
original bytes, hashes and resource ceilings in new directories. The molecular
inspection script performs original-byte/producer/ZIP/XML-shape checks without
creating a processed validation observation. Lane reports document their exact
commands and include unsuccessful retrievals. No source figure was digitized.

No AMPS build, particle run or physical/observational qualification test is
performed. Acquisition and source QA are not I21–I24 PASS. No raw canonical
data, installed database, acceptance/reference, roadmap or Git index is edited.
M1/M2/M3 formal status remains as recorded in STATUS.md.

The source roadmap SHA256 remains
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the data roadmap SHA256 remains
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source resource/development exceptions govern work order. Original
scientific, data-interpretation, calibration/hold-out and predecessor contracts
continue to govern qualification; the two roadmaps are not silently reconciled.
