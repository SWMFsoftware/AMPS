# M3 literature source access

The archive search followed the I22/I23/I24 gaps in the
[r6 research report](../../srcMoon/test/M3_ARCHIVE_SOURCE_SEARCH_20261124.md).
This follow-up uses the existing offline packet importer, shared CSV
Catalogue and versioned installation. The numerical sources are the
published text and table cells in specifically licensed articles.

## Source and test contract

| Source | Numerical view | Scope and missing quantities |
| --- | --- | --- |
| Dibeler, Walker and Rosenstock (1966), DOI 10.6028/jres.070A.038 | Two measured water onsets in eV and Å. Three measured relative ion signals at the paper's 584 Å, 21.23 eV line. | Historical onset convention. Relative signals retain the parent normalization and the original absence of mass/structure response corrections. A fragment mass does not identify its neutral co-product. No plotted spectral values, solar branch probabilities or fragment kinetic partition. |
| Kim (1992), DOI 10.6028/jres.097.032 | Table 1 constants for electron-impact ionization of neutral H and He. | Source fit thresholds and Eq. (2) for this integral cross section. Preserve the coefficient text and the expression convention. The original also contains ionic rows, which are outside the neutral-target projection. No replacement for D05, excited residual-state fraction or covariance/slope bound. |
| Cui et al. (2025), DOI 10.3390/atoms13090074 | Table 1 percentages for capture into ground and excited He, at five projectile energies, with measured and theory cells separate. | Relative single-capture branching for a helium projectile. Keep printed uncertainties and the blank excited measurement at 0.5 keV/u. Raw, absolute angular data and unmeasured tails are still needed for angular qualification. |

The two NBS/NIST Journal of Research articles are individually pinned as US
government technical publications. Keep article citation and requested NIST
acknowledgment. The original Cui article and Table 1 are tied to its printed
CC BY 4.0 notice. These decisions do not cover a different publisher, a
third-party array from a cited archive, or NIST SRD generally.

The recipe reads the native PDF text layer with `pdftotext`, addresses the
specified pages, and retains the published numeric tokens. It refuses changed
native hashes and unexpected tables, identities, cells or source prose.
No plot digitization, numerical fitting, extrapolation or unit conversion is
performed. Source bytes and the exact recipe are independently replayable
from the snapshot and acquisition manifests.

The oracle for an accepted catalogue view is its original source table or
passage, verified independently against the installed API. Test an
article-specific grant's rejection of mismatched source identities in
the shared software suite. On authentic data, compare the published tokens,
source roles, missing cell, row counts, and C++/Python queries. A missing
optional original is SKIPPED, a mismatch is FAIL, and infrastructure or
malformed data is ERROR. These are source-access checks. I22–I24 production
activation and their formal data, numerical and predecessor criteria remain
outside this work package.

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The external roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source development/resource amendment governs development order. Shared
science criteria, data interpretation, immutable originals and predecessor
gates control qualification.

## Archive and publication findings

- **I22 water products.** The original NBS paper gives H₂O⁺ at 12.593 eV
  and the OH⁺ onset at 18.05 eV. Its 21.23 eV observations give relative
  H₂O⁺, OH⁺ and O⁺ ion signals of 1.0, 0.4 and 0.05. These are historical
  uncorrected mass signals. Bodi et al.'s threshold-electron coincidence
  measurement gives the separate 0 K OH⁺ appearance energy
  18.1182 ± 0.0009 eV. Retain the difference in measurement convention.
  The search recovered no complete single-photon electron and fragment
  distributions for H₂O04–06. The I22 solar state crosswalk and incident-OH
  sticking and accommodation are still required. The DOE-hosted
  Samson–Starace report documents He⁺(2s,2p) fluorescence and pressure
  quenching, without a usable broadband state-yield table for D04.
  [Dibeler et al.](https://doi.org/10.6028/jres.070A.038)
  [Bodi et al.](https://doi.org/10.1039/C4SC00791C)
  [Samson and Starace](https://www.osti.gov/servlets/purl/5545053)
- **I23 helium products.** Kim's Table 1 is now available as an alternative
  integral reference. The He fit uses the printed 24.58 eV threshold.
  Bellm et al. measure differential ratios by residual He⁺ principal state
  for particular electron geometries and incident energies, rather than
  complete state- and energy-integrated yield for the 10 eV Maxwellian.
  Renormalizing their ratios as solar or plasma branching would change
  their meaning. NIST SRD 107 includes ejected-electron energy spectra,
  but the NIST general reuse page separates SRD products from the
  technical-publication grant. The retrieved SignalHelium deposit contains
  a mixture of the authors' simulation and embedded third-party
  cross-section files. That mixture remains metadata-only. The
  source-to-state and approximation-error bounds are still required.
  [Kim](https://doi.org/10.6028/jres.097.032)
  [Bellm et al.](https://doi.org/10.1103/PhysRevA.75.042704)
  [SRD 107](https://doi.org/10.18434/T4KK5C)
  [SignalHelium](https://doi.org/10.5281/zenodo.14173434)
- **I24 charge exchange.** Cui's ground and excited branches give new
  comparative evidence for He⁺–He single capture. At 5 keV/u the
  measured excited fraction is 3.11 ± 0.27 percentage points. The missing
  0.5 keV/u excited measurement stays blank. The public angular figures
  use peak normalization. The authors identify an on-request route to raw
  data. The independent Wang CDS catalogue was located. Native table
  downloads timed out or failed, so its numeric cross sections and rates
  are absent from this release. The Groningen institutional catalogue
  lists the sodium thesis and chapters. No native sodium ground-channel
  angular bins or confirmed D119990 attribution were recovered in this
  pass.
  [Cui et al.](https://doi.org/10.3390/atoms13090074)
  [Wang et al.](https://doi.org/10.1051/0004-6361/202244598)
  [Groningen repository](https://research.rug.nl/en/publications/electron-dynamics-in-ion-atom-interactions/)

## Acquisition and replay

The research artifacts are under
`test_output/chemical_database/m3-source-r7-xBR6k1oF/`. Top-level and
`followup/` acquisition plans and receipts preserve native URLs, timestamps,
response status and SHA-256. Failed requests and the original analysis
remain in the research directory. Four held entries retain metadata
and original byte identities. No array deserialization from the
SignalHelium ZIP was needed.

```bash
python3 -B utility/chemical_database/native_m3_literature.py \
  --acquisition test_output/chemical_database/m3-source-r7-xBR6k1oF \
  --output-dir test_output/chemical_database/NEW_M3_LITERATURE_PACKET
```

The four numerical views contain 12 published rows. The native
SHA-256 of the NBS water PDF is
`caab1a72ff76e4d02fe59c559efd2eff4f2a5542730505eb7fc2776dcc9eff04`,
the Kim PDF is
`b777277069851715e5810c17f36ad8fcbaf005fd6c0071c2b9baba28082be961`,
and the Cui PDF is
`7ed7cbf026341d5c1c7be67c25a48000f0b8243a366ba9c5d259040badb922d0`.
`build-command.json` records the complete r6 replay and one additional
r7 packet.
