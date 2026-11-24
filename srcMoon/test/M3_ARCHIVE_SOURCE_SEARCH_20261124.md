# M3 physics archive search — 2026-11-24 UTC

The search confirmed an open numerical angular archive for resonant helium
charge exchange, verified the live IAEA sodium partial-capture table and
located the strongest same-reaction sodium angular measurements. It also
recovered input tables for bounded helium photoionization and electron-impact
studies. The current milestone still needs the charged-water product laws,
the photon-state crosswalk and the electron-impact approximation-error
envelope. The [M3 completion plan](M3_COMPLETION_FOLLOWUP_20261124.md) and
[sodium source review](integration-tests/I24/SOURCE_ADMISSION_20261124.md)
define the cases these inputs could support.

## I24 charge exchange

### Proton capture from ground-state sodium

- The live IAEA D119990 page contains the same 12 numerical knots at
  4–25 keV/u as the r5 source. I reproduced the byte-identical numerical
  projection using the existing decoder's renderer, shape, units and frame
  checks. The live HTML has a different hash and the native download route
  returned HTTP 500. The live page labels H(1s) and Na⁺(2p⁶), gives a
  target-rest projectile energy axis, and identifies a JAEA/Kyoto
  compilation. The original experimental citation, target preparation,
  uncertainty and angles are still required for admission. The prior
  frozen HTML and projection remain the references.
  [IAEA D119990](https://db-amdis.org/collisiondb/datasets/119990/)
- Knoop's Groningen MOTRIMS thesis, chapter 5, explicitly studies
  H⁺ + Na(3s) over 4–25 keV/u. Table 5.1 and the recoil spectra separate
  the outer-shell H(1s) channel from inner-shell capture that also produces
  H(1s) and leaves excited Na⁺. Section 5.2.3 gives outer-shell ground-state
  differential extractions from 7 keV/u upward. Figure 5.9 is peak-normalized
  transverse momentum and Figure 5.10 is peak-normalized relative angle.
  The thesis quotes a 10–20% normalization systematic for the integrated
  partials in Figure 5.6. This systematic alone cannot bound the angular
  tails. The match of reaction label and overall energy range makes
  these measurements a candidate D119990 crosswalk. Establishing the
  lineage requires the original table and a producer or catalogue citation.
  [Groningen thesis, chapter 5](https://pure.rug.nl/ws/portalfiles/portal/2855778/c5.pdf)
  [Zapukhlyak et al., experiment and calculation](https://doi.org/10.1088/0953-4075/38/14/003)
- The full Courbin paper is now available locally. Figures 7–8 give
  center-of-mass angular calculations at about 0.5 and 2 keV for
  capture into H(n = 2), including the initial Na(3s) case. Richter et al.
  measure the H(n = 2) density matrix at small angles. These are appropriate
  leads for an excited-hydrogen case including its state and radiation
  account. Nagata's 0.4–5 keV measurements give effective eventual
  ground-state production inferred from a separate metastable fraction.
  Prompt H(1s), eventual H(1s) and H(n = 2) require distinct interpretations.
  [Courbin et al., archived article](https://oa.upm.es/61369/)
  [Richter et al.](https://doi.org/10.1088/0953-4075/26/4/012)
  [Nagata](https://doi.org/10.1143/JPSJ.48.2068)

For the approved sodium case, request the native outer-shell H(1s) energy
and differential bins from the Groningen authors or the IAEA data
contributors. Specify the Na(3s) preparation, Na⁺(2p⁶) recoil selection,
momentum-to-angle conversion, acceptance and efficiency, uncertainties,
normalization and coverage of the tails. If only a bounded-error forward
law is available, it still needs a source-supported angular tail bound and
the existing heavy-product velocity and energy error review. A plotted
angular extent supplies no tail bound.

### Resonant helium as a separate controlled case

Swaczyna et al.'s CC BY 4.0 Zenodo deposit contains a separate
`He0He1ct` numerical branch. The native size and MD5 matched and the ZIP
integrity check succeeded. I parsed all 661 energies from 0.0001 to 10,000
eV in the center-of-mass frame, each with 768 angular rows. The second
column already includes `2π sin(θ)`. Summing that column times the
provided quadrature weights reproduced every printed integrated cross
section with maximum relative discrepancy `4.60e-11` (rounded up).
At 10 eV center-of-mass, the integral is `2.081618e-15 cm²`.
Appendix A defines the neutral scattering angle for which the exchange
term peaks near π. A mover using the projectile direction must apply the
corresponding product-identity and angle transformation.
[Zenodo numerical archive](https://doi.org/10.5281/zenodo.4555716)
[Swaczyna et al., Appendix A](https://pmc.ncbi.nlm.nih.gov/articles/PMC8861972/)

This is a candidate for a controlled He⁺ + He ground-channel study, subject
to review of the identical-particle treatment, state and energy window,
and comparison to independent rates and angular evidence. Wang et al.
provide an independent state-selective theory reference in eV/u. The CDS
catalogue endpoints timed out during retrieval. The integrated angular
check is a source-format check. It supplies neither an independent
physical error estimate nor a result for proton–helium or proton–sodium.
[Wang et al.](https://doi.org/10.1051/0004-6361/202244598)
[CDS catalogue](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/674/A201)

## I23 electron-impact product model

Curtin's CCC archive supplies 55 ground-parent He single-ionization
rows in `TICS.s1S`, covering 24.5–2000 eV incident electron energy.
The NIFS evaluated report supplies a separate fitted cross section in
Eq. (9) and Table IV, based on experimental and CCC data. Rejoub et al.
measured charge-resolved He, Ne and Ar ionization from threshold to
1000 eV. Charge-state yields leave the populated internal states
unresolved. The NIST SRD 107/author BED energy distributions and
laboratory measurements of differential electron ejection at 100 and
200 eV can constrain restricted partition checks.
[Curtin CCC archive](https://atom.curtin.edu.au/ccc-db/index.html)
[NIFS assessment](https://www.nifs.ac.jp/report/NIFS-DATA-059.pdf)
[Rejoub et al.](https://doi.org/10.1103/PhysRevA.65.042713)
[NIST SRD 107](https://www.nist.gov/pml/electron-impact-cross-sections-ionization-and-excitation-database)
[Differential laboratory measurements](https://doi.org/10.1103/PhysRevA.16.26)

For the existing stationary ⁴He, 10 eV Maxwellian pilot, I integrated
`E σ(E) exp(−E/Te)` using linear interpolation strictly within the
tabulated Curtin domain. Halving the quadrature step bound from 0.1 to
0.05 eV changed the incident mean by about `1.3e-9 eV`. The
finite-table diagnostic gives `7.21074e-16 m³/s` and `46.08014 eV`
conditional incident mean. The frozen Voronov values are
`7.76063e-16 m³/s` and `45.41818 eV`. These differences are
diagnostics of source/model spread. About 9.81% of this diagnostic
ionization integral lies above the first He⁺ excitation onset
(approximately 65.4004 eV incident energy). This percentage measures
energy eligibility in the tabulated model. The actual excited-He⁺ yield
and its uncertainty remain unresolved.
[NIST helium ionization energy](https://www.nist.gov/system/files/documents/srd/jpcrd690.pdf)
[NIST He II energy levels](https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable6.htm)

The pilot still needs evidence for the energy dependence and uncertainty of
the cross section or Maxwellian rate slope, tail contributions and
state-resolved products or an independently bounded excited fraction. A
shared absolute normalization error cancels from the moment ratio. It
cannot constrain its energy-dependent error. Density and temperature
prescriptions remain verification inputs, with measured D06 forcing
required for corresponding physical claims. The public CCC README labels
the current tables unpublished and the earlier rights review found no
managed-storage grant.
[Source rights review](../../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261119.md)
[Existing model-error contract](M3_RESERVOIR_PRODUCTS_CONTRACT.md)

## I22 photon products and U24 surface response

- **Atomic products.** The NORAD He ground-parent entry contains 1,560
  native points for ground residual He⁺ and 1,777 in the separate total
  table. There are 1,467 common energies and 208 common energies with
  unequal cross sections. The native zero-cross-section threshold entry
  is 1.773760 Ry in both. It is a source energy, with a different value
  from the spectroscopic threshold. The recovered Verner–Yakovlev data
  identify Ne and Ar subshell ionization. Their shell labels do not
  identify individual residual-ion levels or the satellite population.
  A finite photon band below the first alternative product-state
  threshold can isolate a ground-state event check. Broadband D04 still
  needs a documented common energy axis, spectrum weighting and residual
  states or a bounded reduced law at both solar endpoints.
  [NORAD ground-residual table](https://norad.astronomy.osu.edu/he1/he1.ptpx.ls.txt)
  [NORAD ground-parent totals](https://norad.astronomy.osu.edu/he1/he1.px.gd.txt)
  [Verner–Yakovlev native tables](https://www.pa.uky.edu/~verner/photo.html)
  [PHIDRATES](https://phidrates.space.swri.edu/)
- **Charged water.** The NIST assessment collects OH⁺, O⁺ and H⁺
  marginal yields. Bodi et al.'s threshold-electron coincidence
  measurement gives `18.1182 ± 0.0009 eV` for OH⁺ + H + e⁻ and
  provides a near-threshold breakdown check. The Caltech archive
  contains rotational distributions for near-threshold surviving H₂O⁺.
  The search found no complete single-photon distribution over the
  D04 solar domain for `H₂O04–06`. The outstanding measurements are
  the coincident electron energy, heavy-pair kinetic-energy release and
  angular or momentum distribution, with neutral co-fragment and
  product internal-state assignments at documented photon energies.
  An O⁺ marginal alone leaves its neutral co-fragment unidentified.
  [NIST photoionization assessment, figures 49–54](https://www.nist.gov/system/files/documents/srd/jpcrd339.pdf)
  [Bodi et al.](https://doi.org/10.1039/C4SC00791C)
  [Caltech archive](https://authors.library.caltech.edu/records/0av93-f0g90)
- **U24 and ion–surface response.** The existing gold-relative H₂O
  sticking data for bare olivine and molecular release references
  remain useful within the documented experiments. Published OH
  diffusion on UV-irradiated amorphous water ice measures the behavior
  of OH already on the surface. The search did not recover absolute
  incident-OH sticking and post-impact accommodation for the required
  substrate and temperature range. The Chang’e-6 NILS archive provides
  negative-particle counts summed over surface-pointing directions under
  CC BY-NC-SA 4.0. The associated NILS analysis and the ARTEMIS
  scattered-proton study constrain hydrogen branches under specific
  instrument and incidence conditions. Na⁺ and He⁺
  neutralization, reflection and retention require their own outgoing
  charge, energy and angle measurements and a balance relative to the
  incident flux on the specified surface.
  [Laffon data](https://doi.org/10.6084/m9.figshare.13274840.v1)
  [OH diffusion experiment](https://arxiv.org/abs/2210.13753)
  [NILS public data](https://doi.org/10.5878/v4sy-3t95)
  [NILS observation paper](https://www.nature.com/articles/s43247-025-02399-7)
  [ARTEMIS proton scattering](https://doi.org/10.1029/2018JA025486)

## Provenance and verification scope

Source revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a`.
The source [roadmap](../development_and_validation.md) hash is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The external roadmap at `/data/vtenishe/moon_validation_data` hashes to
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source tree supplies the early development and resource amendments.
The canonical shared requirements govern formal science criteria,
preserved data and the I13–I20 predecessor gates. The exact r5
`database-lock.json` hash is
`7b588795ba3dcbd17491d5afbb9a417f62824970afa94b10645e8eda835799a2`.
D04, D05 and r5 references were used in their existing scope.

The local evidence directory is
`test_output/srcMoon/m3-source-search-9wqcZyI0/`. The root
`retrieval.json` and the `followup/`, `additional/` and `confirmed/`
manifests record URL, byte count, SHA-256 and request outcome. A successful
HTTP request only authenticates the received bytes. The Nagata nominal PDF
route returned HTML, the attempted PMC mirror resolved to a different
article, and the Bodi supplement route failed in the shell. The report
cites the verified publisher records for those papers.

| Verified artifact | SHA-256 |
| --- | --- |
| Live D119990 HTML | `1015ec7a843e02ef7fc1ddf97f3589ac4dc4823f12bb28e9434313876dd90fd6` |
| D119990 numerical projection, identical to r5 | `2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93` |
| Groningen thesis chapter | `7daca10c4ecf37bd79f0a5d353cfc9c6d6adc9c2b0277d74aa01d5381c5c1785` |
| Zenodo numerical ZIP | `16027a414c603d84a2a500964a8a877a2bf2e8d5573f53827330f75e5ab7e0e8` |
| Curtin He numerical ZIP | `e4a3564c1d6f811134761496cf00c3e7026d26fbb4d6d9748c72230e34e8ecb3` |
| NORAD He ground residual | `d89955666852f3394dd7edaf8f5cd65b86d8b2d9d8a853976b6b45973a698615` |
| NORAD He total | `d51f995f3beed2367cebf25e8e90e4e240e16e685a429609bb6a1c28542f5074` |

From the AMPS root, the research-only audits can be reproduced with
the following commands against the saved bytes.

```bash
evidence=test_output/srcMoon/m3-source-search-9wqcZyI0
unzip -tq "$evidence/confirmed/native/zenodo-he-diff-cross-sections.zip"
python3 -B "$evidence/audit_native.py" "$evidence"
python3 -B "$evidence/audit_zenodo.py" "$evidence"
```

This package ran zero AMPS builds or linked tests. The last formal
I21–I24 outcome remains 0 PASS, 0 FAIL, 0 ERROR and 4 SKIPPED. Database
admission and any changed I24 channel require separate review under the
existing contracts. The I21 constant-hazard numerical check and the rest
of its qualification matrix remain open independently of these data.
