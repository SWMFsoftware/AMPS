# Chemical database search and Hartley2 H2O reuse audit

Research date: 2026-11-18. AMPS revision:
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`, including the previously recorded
uncommitted Moon development. This package inspects sources and downloads
research evidence. It changes no production physics, input selection, accepted
reference, numerical criterion or canonical validation data. Formal M3 remains
**SKIPPED / NOT QUALIFIED**.

The practical finding is that suitable **rates are available**, and a reachable
H2O/OH network and recoil implementation already exist in the comet model.
Reuse still needs an explicit lunar product/state/energy contract and correct
physical-time integration. A total rate, mean energy and conservative numerical
sampler establish different pieces of that contract.

## Authority and retained evidence

Source roadmap SHA-256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
External roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source includes approved early-development/resource amendments; common
formal U/I meanings, predecessors and thresholds remain binding. Native data
metadata control scientific interpretation. Neither roadmap was edited.

Evidence directories relative to the repository root:

- `test_output/srcMoon/m3-database-hartley-research-20261118T083645Z/`:
  consolidated source audit, native-paper cross-check, baseline/protection
  identities and documentation correction.
- `test_output/srcMoon/m3-hartley-h2o-code-audit-20261118T083318Z/`:
  `audit.md`/`audit.json`, exact source snippets, production paths and hashes.
- `test_output/srcMoon/m3-molecular-databases-20261118T083551Z/`:
  `REPORT.md`, acquisition manifests, native molecular files and
  `HDF5_INSPECTION.md` with its reproducible read-only inspection.
- `test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z/`:
  `REPORT.md`, `native_inventory.json`, `packet_manifest.json`, source queries,
  numerical extraction scripts and exact successful/failed retrieval records.

All fresh downloads remain research evidence. They were not promoted into a
canonical Dxx package, redistributed with production code or used for a
validation PASS. Native responses, including unsuccessful ones, are retained.

## Existing solar baseline: PHIDRATES

[PHIDRATES](https://phidrates.space.swri.edu/) provides separate solar cases,
channel identities, integrated rates and average excess energies. The frozen
local D04 package contains 15 reactive channels: four atomic, seven H2O and
four experimental OH channels. Its two absorption-total records are not extra
destructive events. Live CGI access timed out during this search; the existing
immutable responses remain available.

The following are actual archived values at **1 AU, optical depth zero**:

| Channel | Quiet rate, s^-1 | Active rate, s^-1 | Quiet / active mean excess, eV |
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

The reproducible audit matches all **30 rate/mean-energy records** to the
frozen physical-specification JSON and verifies all **11 processed-product
hashes**. `frozen_phidrates_rates.csv` in the consolidated evidence includes
the atomic rows and exact source locators. No uncertainty or product-state
population was inferred from the printed precision.

The downloaded author-hosted
[Carberry Mogan et al. 2022 paper](https://wurz.space.unibe.ch/Mogan_JGR_2022.pdf),
Table B1, supplies the seven active-Sun water channels at 5.2 AU. Reading its
native PDF text and applying its stated inverse-square normalization produces
rate intervals compatible with all seven archived active values and the same
mean energies. This checks transcription and normalization; both compilations
derive from Huebner, so it is not an independent physical benchmark. The paper
lists some O products without resolved state; numerical matching does not
qualify their state mapping. Its electron-impact rates describe a particular
Callisto electron population and are not universal lunar coefficients.

One existing documentation transcription error was corrected: the physical
specification Markdown printed Ar active as `9.07e-7 s^-1`; its frozen JSON,
native response and CSV consistently contain **7.74e-7 s^-1**. Only that Markdown
cell changed. The previous document and before/after hashes are retained in
`documentation-correction/`; no production rate or data file was changed.

## Solar-activity alternative already used by I3

[Moore and Mapaye 2025](https://doi.org/10.3847/PSJ/add6a2) supplies
`j_1AU = A0 * F10.7P^alpha`, with the proxy combining daily and 81-day solar
flux. Its fits use 2004–2014 irradiance. Two examples are H2O→H+OH
`A0=3.60e-6, alpha=0.24`, and H2O→H2O++e
`A0=4.23e-8, alpha=0.56`. These rows were checked against the indexed primary
table. Publisher retrieval returned a CAPTCHA, not a usable article download;
a complete native table/uncertainty cross-check remains necessary for import.

The comet's existing `srcHartley2/i3/reference/photochemistry_channels.csv`
contains 17 selected coefficients; `prepare_i3_cases.py` generates positive
per-case values in `input/i3.in`. These dates, solar indices and comet source
conditions must not be copied as lunar observations. A Moon implementation
would declare its own radiation model, proxy input and distance normalization.
The source-law evaluation is reusable; lunar runtime wiring is separate work.

## Other online sources: concrete use and limits

| Database / source | What was acquired | Appropriate role |
|---|---|---|
| [Leiden photochemistry](https://home.strw.leidenuniv.nl/~ewine/photo/) | Rate CSV, quiet-Sun spectrum, H2O/OH cross sections and water partial HDF5 | Spectral and aggregate-branch cross-check; source coverage/state limitations remain explicit |
| [MPI-Mainz atlas](https://www.uv-vis-spectral-atlas-mainz.org/uvvis/) | Published water absorption datasets with temperature and wavelength metadata | Absorption references; require quantum yields and a radiation field for reaction rates |
| [ExoMol](https://www.exomol.com/data/data-types/photo/) | OH MYTHOS and H2O DTU numerical spectra/metadata | Temperature-dependent absorption; not complete outgoing particle distributions |
| [UMIST](https://umistdatabase.uk/downloads), [KIDA](https://kida.astrochem-tools.org/networks) | Native reaction-network files and selected photon rows | Network audit; their ISRF rates and reduced products cannot be treated as solar, state-resolved channels |
| [NIST electron-impact database](https://physics.nist.gov/PhysRefData/Ionization/) | H2O BEB totals, orbital data and He singly differential energy tables | Molecular electron total/energy-sharing references, with declared beam-energy support |
| [OPEN-ADAS](https://open.adas.ac.uk/) | Atomic Maxwellian coefficients, metastable He rates and He capture tables | State-labelled rate alternatives and partial capture data; preserve source assessment and format conventions |
| [LXCat / original Phelps source](https://nl.lxcat.net/notes/index.php?download=phelps5) | Original Na3s/3p ionization tables, plus selector metadata | Beam total cross sections; restricted source domains and channel references remain binding |
| [IAEA CollisionDB](https://db-amdis.org/collisiondb/) | Native Na partial-capture measurements and separate composite-fit metadata | State-specific cross-section reference, without inventing angular or all-state closure |

Leiden's [radiation-field documentation](https://home.strw.leidenuniv.nl/~ewine/photo/radiation_fields.html)
specifies that its normalized Solar rates require **×37700** for approximate
1-AU solar intensity. The corresponding dissociation values are H2O
`8.671e-6 s^-1` and OH `7.2384e-6 s^-1`. Preserve the original normalized
values separately. This is a documented conversion, not a calibration factor;
the spectrum and wavelength coverage differ from D04. Raw ISRF coefficients
must not receive an arbitrary solar multiplier. The retained spectrum covers
67.057–419.975 nm and is not complete short-wavelength EUV support.

The native Leiden water HDF5 was inspected in an isolated `/tmp` Python
environment. It has seven arrays with 21,723 samples, but no resolved
daughter-state/energy/angle distributions. It also has 343 negative samples
in each OH/O/total-dissociation array and nonclosing partial sums. At the
largest absolute partial/total residual, the discrepancy is 3.41% of the
local total. No clipping, mask invention, splicing or renormalization was
performed. Its use as a production hazard table needs source interpretation.
ExoMol also has a native metadata/table row-count discrepancy recorded in the
molecular report. These findings do not invalidate the separate frozen D04
integrated-rate baseline.

The NIST H2O BEB table covers 12.61–5000 eV; at 100 eV its total cross section
is 2.348 Angstrom². The He energy-differential tables do not supply correlated
electron angles. Beam cross sections need an explicitly specified distribution
before forming an electron rate. No unsupported tail extension is selected.

OPEN-ADAS has an evaluated He capture file spanning 0.3–500 keV/u, including
state information at ordinary-keV energies. Its low-energy fractions and
high-state completion include source model assessments, and joint angles are
absent. The whole-element Dere files also supply neutral Na/He/Ne/Ar Maxwellian
ionization rates; a charge-zero-only web query misses them. These alternatives
do not replace the currently selected Voronov table automatically. Native ADAS
redistribution terms need review before packaging data with distributed code;
the current copies are local research evidence.

The original Phelps compilation provides 16 Na3s ionization nodes through
100 eV, threshold 5.14 eV, and a separate 12-node Na3p table. Current numerical
export access was challenged; the separately acquired original tables are the
actual numerical source. IAEA D119990 contains 12 measured H(1s)+Na+ capture
points at 4–25 keV/u. **D107581 combines capture and proton-impact ionization**;
it is not a replacement for a pure capture total. Native embedded numeric
arrays were extracted directly, with fitted lines and plot-band polygons
kept distinct; no curve digitization was used.

## Exact Hartley2 reuse boundary

The canonical `input/hartley2.input:153` sets `ChemistryMode=off`. The reachable
I3 configuration instead selects seven neutrals and 17 explicit channels:

```text
SpeciesList=H,H2,O,OH,H2O,CO,CO2
ChemistryReaction=H2O_H_OH | H2O | H + OH | HELIOCENTRIC | 0 | 0 | 0 | 1 | 4.04
ChemistryMode=on
```

The zero coefficient is a compile-time placeholder. Required runtime records
replace it before allocation. An empty source-tree generated header is a
sentinel, not evidence that the configured network has no water chemistry.
`CometRuntimeInput.cpp`, `CometChemistry.cpp` and the shared
`ChemicalSolver::cReactionNetwork` are the actual production path.

There are two distinct comet model choices:

| Model | Three neutral water rates at 1 AU | Declared excess energies |
|---|---|---|
| Rubin2011 catalog in `input/CometPhysicsDatabase.pm` | `1.00e-5`, `6.00e-7`, `7.50e-7 s^-1` | `1.70`, `3.80`, `0.70 eV` |
| I3 runtime network | Case-dependent activity law | Fixed compiled `4.04`, `3.94`, `0.697 eV` |

The Rubin catalog's additional ionization row is `NONE`, and I3 has nine
neutral-system sink rows. These remove neutral particles; they do not create
the physical ions/electrons required by Moon M3. State-labelled O channels
also collapse into one compiled O species in this comet configuration.

`PhotolyticReactions::Huebner1992ASS::GetProductVelocity_2` and the corresponding
comet two-body construction provide useful opposite-momentum, isotropic recoil
math at a declared translational energy. This is a candidate **reduced physical
model**, subject to review of energy meaning, internal states, isotropy,
photon momentum and validity. A measured full joint PDF is not the only possible
model contract; an approximation must be explicitly supported and bounded.

The many-body comet sampler draws Gaussian velocities, removes their mass-weighted
center-of-mass motion and rescales the energy. It is not generally uniform
microcanonical sampling for unequal masses. Its energy partition needs a
separate model decision. I3 changes rates with solar activity but retains
fixed recoil energies; its water energies correspond to the active table above.
This should not become an implicit quiet-Sun energy choice.

Do not copy these comet executor details into Moon unchanged:

- Stationary `W/dt` carriers and daughter timestep ratios: Moon uses represented
  physical-number packets and retained particle intervals.
- Whole-step loss sampling: Moon needs its actual event/exposure chronology.
- Parent auxiliary/internal-state cloning: daughter states need explicit treatment.
- Legacy approximately 0.2% mass-table tolerance: model mass and formation-energy
  conventions must be consistent, not silently loosened for binding energy.
- `NONE` ionization products: Moon needs complete tracked/bath conservation.

The code audit records exact lines/snippets and existing parser, core and
ensemble-test entry points. Those tests were not run during this research.
The comet solver's DIM1/domain/source and collision configuration is not a
validated lunar model or an independent observation.

## Recommended next work, requiring model review before activation

1. Preserve the current D04 quiet/active source as the explicit baseline.
   Reuse the available rates; no new invented rate is needed for these channels.
2. Consider the existing I3 solar-proxy law as a separately named optional model,
   with complete primary coefficient/uncertainty verification and lunar inputs.
3. Prepare a concrete two-body neutral-water closure from the existing recoil
   math and a consistent chosen source energy. Review its state, energy and
   angular approximations before implementing physical particle production.
4. Use Moon's conservative controller, physical-time packets and force tracers;
   reuse comet mathematics/catalog interfaces without importing steady-flux
   resampling or neutral-only sinks.
5. Qualify new electron/capture candidates independently and retain missing
   angular, state, surface and source-domain coverage. Three-body/cascade and
   ion processes need their own complete models.

No physical model choice was approved or activated by this search. This
research narrows the input gaps; it closes no U/I gate or formal milestone.

The subsequent [concrete decision proposal](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md)
uses these findings to specify P01–P09 recommendations, explicit reduced-model
approximations, an ideal verification collector and a source-specific He
laboratory benchmark. All decisions remain pending. That proposal does not
change this research's data-readiness, source identities or qualification scope.
