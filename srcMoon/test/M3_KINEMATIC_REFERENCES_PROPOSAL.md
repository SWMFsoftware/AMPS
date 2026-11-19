# M3 daughter-kinematics sources and acquisition proposal

Research update, 2026-11-17. **Source suggestions only: no kinetic model is
approved or activated by this document.** No observation, coefficient, frozen
test threshold, channel membership, canonical data file or physical default
was changed. The approved conservative event executor can accept qualified
event products; it cannot obtain missing distributions from a total rate.

## Recommended acquisition order

1. Obtain channel- and photon-energy-resolved Na photoionization cross
   sections **and angular information**, together with the solar photon
   spectrum and its exact averaging convention. Start with the Na-specific
   sources below, including their Cooper-minimum limitation. Keep the
   existing D04 total-rate checks separate from this kinetic qualification.
2. Obtain the Na ground-state electron-ionization joint differential data
   from the Na CCC calculation, covering the incident energies needed by the
   declared Maxwellian bath. A few fixed scattering geometries or a total
   cross section cannot determine the two outgoing-electron distribution.
3. For molecules, first obtain the wavelength/state-resolved H2O→H+OH
   distributions and the initial H2O-state contract. Then obtain the other
   complete D04 channel distributions individually. Retain the full original
   rate family even when a kinetic channel is unavailable; refuse requested
   transport rather than delete a hazard or redistribute its probability.
4. Review the acquired numeric files and their units, frames, energy zeros,
   validity ranges, interpolation and polarization conventions before
   selecting an executable kinetic policy. No author was contacted in this
   work. Contacting authors remains a separate, explicit action.

## Primary sources and what each can establish

The distinctions in the final column are requirements of the proposed AMPS
application, not a claim that the papers attempted to provide every needed
quantity. Author-request availability is an acquisition route, not data
already present in this checkout.

| Process | Primary source | Concrete acquisition route and limitation |
|---|---|---|
| Na valence photoionization | [Theodosiou2022, arXiv2212.00577v1](https://arxiv.org/abs/2212.00577v1) | Open versioned PDF; ground/excited cross sections around threshold/Cooper minimum. This is useful rate/spectral evidence, not by itself a full solar-weighted electron-angle sampler. Ground-state measurement/theory disagreements are discussed. |
| Na photon/atom calculations | [Singor, Fursa, Bray & McEachran2021, Atoms9,42](https://doi.org/10.3390/atoms9030042) | Open article and stated author-request numerical data. Relativistic/quasi-one-electron calculations provide a second Na-specific cross-section source. Ask for numerical channel amplitudes/angular conventions, not only figure coordinates. |
| Na photoelectron angles | [Hosea etal2023, Atoms11,125](https://doi.org/10.3390/atoms11100125), [open OSTI PDF](https://www.osti.gov/servlets/purl/2345057) | Explicit Na3s dipole/quadrupole angular formalism and Cooper-minimum analysis; numeric data are author-request. It demonstrates why a universal β=2 or isotropic prescription should not silently be applied over all photon energies. |
| Na electron-ionization joint distribution | [Bray, Fursa & Stelbovics2008, J.Phys.B41,215203](https://doi.org/10.1088/0953-4075/41/21/215203) | Na-specific fully differential CCC calculation. Request the underlying incident-energy/outgoing-electron grids or amplitudes, spin/normalization rules and uncertainty/convergence evidence. Full publisher content/numeric data were not acquired here; do not claim its coverage fills the required Maxwellian support yet. |
| Na electron-impact numeric access | [Curtin atomic CCC database](https://atom.curtin.edu.au/CCC-WWW/) | Official selector exposes NaI ground3s and total ionization (`TICS`), with cm² as an explicit option. The selector was retrieved; its final TICS CGI request returned HTTP500. No numerical Na table was recovered. Its total cross section alone would still not supply recoil. The database also warns that inner-core ionization may matter. |
| H2O→H+OH state distribution | [Tabone etal2021, A&A650,A192](https://doi.org/10.1051/0004-6361/202039549), [open versioned text](https://arxiv.org/html/2101.01989v2) | §2.2/AppendixA collect electronic/vibrational/rotational OH probabilities versus wavelength from vanHarrevelt/vanHemert calculations, including repeated calculations with an updated surface. Request these numerical probabilities and original calculation metadata. They concern rotational-ground-state H2O, omit OH contributions below114nm and contain stated channel approximations; they are not a complete solar H2O kinetic law. |
| H2O product-energy/angle evidence | [Chang etal2019, NatureCommunications10,1250](https://doi.org/10.1038/s41467-019-09176-z) | Open PDF and supplementary information; H-atom time-of-flight/polarization measurements at121.57,117.5,115.2nm provide state/energy and angular constraints for jet-cooled H2O. Underlying study data are author-request. Three wavelengths and that initial population cannot silently become a broadband solar/warm-parent joint distribution. |
| OH excited-state dissociation | [vanDishoeck & Dalgarno1983, J.Chem.Phys.79](https://doi.org/10.1063/1.445864), [NASA record](https://ntrs.nasa.gov/citations/19830057250) | Primary state/transition-moment and dissociation-mechanism reference. NASA exposes metadata but no downloadable article for this record. An archival article or author numerical data is still needed before using its quantities. |
| OH modern machine-readable candidates | [Mitev etal2025, MNRAS539,3732](https://doi.org/10.1093/mnras/staf729), [official MYTHOS files](https://www.exomol.com/data/molecules/OH/16O-1H/MYTHOS/) | Open methods and downloadable state/model/continuum files; the publication covers temperature-dependent absorption/pre-dissociation, and state-resolved forms. These are promising independent spectral/state data. A temperature-averaged absorption table is not automatically a state-resolved O+H angle law, and substituting it for frozen D04 `OH_exp` needs explicit review. |

## Why rates do not close the kinetic problem

These are physical requirements inferred from conservation and conditional
sampling, not numerical values taken from one paper.

For photoionization, select a photon energy/channel according to the
declared incident spectral photon flux times the corresponding partial
cross section. Conditional outgoing-electron energy, direction and residual
ion state then determine ion recoil through total energy/momentum closure,
including the incident photon. The D04 scalar mean excess energy does not
determine this conditional distribution. A dipole/unpolarized approximation
must name its photon-direction axis and scope; it must not be confused with
a laboratory linear-polarization formula.

For electron ionization, event-conditioned incident electrons are weighted
by relative speed and ionization cross section. The two outgoing electrons
and residual-ion state must share the available energy/momentum. The
already selected Maxwellian Voronov coefficient is an integrated rate;
inverting that scalar coefficient cannot uniquely recover either the
collision-conditioned incident energies or the joint outgoing distribution.
Replacing only its kinetic sampler with an unrelated differential source
also requires an explicit consistency/normalization contract. Do not
silently change the approved total rate to hide that difference.

For two-body molecular breakup, specified initial/final internal states and
photon energy determine the relative translational energy and its mass
partition; the direction still requires angular evidence. Three-body breakup
needs joint momentum correlations. Product internal state must remain
explicit because it controls both translation and subsequent chemistry.
Separate `OH(X)`, `OH(A)` and predissociative states cannot silently be
combined into one ground-state kinetic particle. The existing D04 `OH_exp`
identifier describes an experimental compilation, not an initial-state
population or a resolved spectroscopic state.

## Concrete policy decision after acquisition

The recommended first transport policy is **a complete, bounded, explicitly
declared state/spectral domain whose required joint data are available**.
That domain must be chosen for source validity before examining model
results. Outside it, either provide the remaining qualified channels or
reject requested kinetic transport. A rate-only closed-box run can still
test the complete original rates, branches and bookkeeping. It cannot
qualify ion pickup, daughter escape, LOS brightness or molecular migration.

The decisions still needed are the actual initial-state populations,
supported photon/bath domain, numerical datasets, energy/mass convention
and uncertainty/normalization policy. The user supplied no preferred
references; this document supplies candidates rather than asking the user
to invent missing physical values. No fixed recoil approximation is
recommended without a quantified, reviewed error bound.

## Production adapter and supporting verification

The implementation boundary is
[`ReactionExecution.h`](../species/ReactionExecution.h). A caller supplies
explicit policy/evidence/frame/closure bounds and each product occurrence:

```cpp
// Fill from a qualified, channel-resolved physical sampler, not a total rate.
ReactionExecution::KinematicPolicy policy;
// policy.qualified, policy.identity, policy.evidenceIdentity, eventFrame,
// maximumClassicalSpeedMS, bounds, bath and products are explicit.
// MakeEventPlan refuses missing policy; PlanScope::RateOnly cannot PIC-commit.
```

The actual field/API definitions control compilation; this excerpt documents
responsibility rather than selecting a production sampler. Products retain
the parent's physical trace interval and get fresh optical clocks. Their
reference-weight correction preserves represented number without a species
birth-time-step factor. Reserve/stage/commit retries must reuse the same
plan and draws. These services have **56 local analytical/control PASS** in
the retained run below; that is not evidence that any physical kinetic
policy above is qualified or that linked I21–I24 PASS.

```bash
/usr/bin/python3 -B srcMoon/test/test_m3_reaction_execution.py \
  --output-dir "test_output/srcMoon/m3-reaction-execution-$(date -u +%Y%m%dT%H%M%SZ)"
```

## Provenance and reproducibility

Inspection revision: `74974c2b2fc1b4cf09a73c320064a22370c0ca94`, with
concurrent uncommitted M3 work. Source roadmap SHA-256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
external roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source roadmap adds approved development-order/status amendments; the
external data contract still controls native data interpretation. Neither
copy makes total rates into kinetic distributions or changes formal gate
order. See the earlier [molecular inventory](M3_MOLECULAR_KINETICS_PROPOSAL.md)
for frozen D04 identities and channel gaps.

Research-only files are under
`test_output/srcMoon/m3-kinematic-references-20261117/`.
`retrieval_manifest.json` records each exact URL, query, UTC retrieval,
byte count, SHA-256, successful download and failed request. Four open PDFs
were downloaded, along with official CCC selector HTML. `pdftotext -layout`
provided searchable local inspection copies; no figure was digitized and
no numeric science table was constructed. The original restricted attempt
could not resolve DNS; authorized host retrieval succeeded. MDPI direct
HTML returned403 and the CCC TICS CGI500; these are acquisition limitations,
not physics test failures. Browser-visible publisher records supplied the
remaining bibliographic/data-availability checks. Canonical validation data
were not modified.
