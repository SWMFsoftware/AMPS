# Atomic state and molecular-product source follow-up

This work uses AMPS revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a`
and continues the archive and publication checks in
[r6](NATIVE_M3_ARCHIVE_SOURCES.md) and [r7](NATIVE_M3_LITERATURE_SOURCES.md).
Public archive records were consulted before the publication text. The
retrieval receipts and source bytes are in
`test_output/chemical_database/m3-source-r8-tU34BqBU/`.

## Preintegration contract

The AMPS roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The validation-data roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The AMPS copy adds approved early development and local verification
exceptions and records practical resource limits. The shared scientific
and dataset-specific prerequisites still control formal I21–I24.

The application-facing entry points for this work are the shared
Catalogue's read-only table and source APIs. The source-packet importer
must bind each newly copied source to its exact bytes and authoritative
rights evidence. The publisher must carry forward the preceding packet,
species, reaction and model catalogue. No AMPS driver, chemistry
profile, default lock, reaction selector or canonical validation data
requires a change. The output is useful independently of the Moon.

The extraction oracle is the original source cell or printed numeric
token identified below. Verify all workbook cells in scope directly
against the XLSX package, and check the published final ratios and
theoretical cross sections independently of the text-extraction recipe.
Check actual CSV access from C++ and Python, and verify that access to
the embargoed record yields metadata only. These are deterministic
database checks with no random seed, MPI or lunar frame.

## Numerical sources

| Source | Exact scope | Native SHA-256 |
|---|---|---|
| [Kalha et al. Figshare v2](https://doi.org/10.6084/m9.figshare.12389750.v2), original XLSX | 103 elements in the original file. The CSV projects the He, Ne, Na and Ar subshell tables from rows 2–17 and preserves cell references and empty cross-section cells. | `95770ea450655b6a26598dfe6b34a64787a558fdb933c1c13db2329da86570b5` |
| [Official Figshare version record](https://api.figshare.com/v2/articles/12389750/versions/2) | Version 2, file 22860461, CC BY 4.0, native MD5 `fdfb61a3867f2965115c5b3d421d8802`. Credit Kalha et al. and [Yeh and Lindau (1985)](https://doi.org/10.1016/0092-640X(85)90016-6). | `cde84aa0eb1d55dcc8f9adb2d2a6b34adf424e2bac089b2efe8a0f0fc97fdaa8` |
| [Sutton, NASA NTRS 19720018903](https://ntrs.nasa.gov/citations/19720018903), original report | Printed p. 31, Table 2 and final estimates. Gas-phase He at 200 eV incident electron energy. Three residual He⁺ orbital cross-section ratios relative to 4s, with author's estimated relative 1σ errors and instrument-width sensitivity. | `d2996c3204e2bc56fa3fe0448bb56c6947e017a17ce5f18a90d62765082904ba` |
| [Official NASA citation record](https://ntrs.nasa.gov/api/citations/19720018903) | Public distribution and government public use permitted. The record identifies the exact report and no third-party material. Cite J. F. Sutton and NASA Goddard. | `ad3ccbff11faba160af354bc12fa0832156d797b3a6b973bc2fcbef6c5ca4871` |
| [Bhatia, Atoms 10, 9 (2022)](https://doi.org/10.3390/atoms10010009), original PDF with CC BY 4.0 notice | The 15 author-hybrid-theory entries for He in Table 2 and eight author's exchange-approximation entries for He in Table 8. Photon ionization of ground-state He. Preserve outgoing electron momentum `k` in the article's Rydberg convention. Table 2 uses Mb. Table 8 gives dimensionless ratios to the ground residual. | `939c68ea8fafa4c8502a126c76b88d2114290cee1095b1c8b875916eef3f1f80` |

The workbook contains tabulated theoretical cross sections in Mb at
discrete photon energies in eV. Its original explanatory sheet gives
the unit and scope. The tabulated subshells do not enumerate all
residual-ion terms or their branching fractions, nor do the energy
knots establish a solar-weighted rate. Bhatia's table preserves its
missing `k=1.2` knot and carries no assigned uncertainty. The CSV uses
the author's calculation column. The article also prints separately
cited experimental and R-matrix comparison columns.

Bhatia's Table 8 gives the theoretical He⁺ 2s-to-1s ratio and
the combined 2s-plus-2p to 1s ratio. Its R2 decimal exponent is
`−1` at `k=0.3` and `0.4` and `−2` at the neighboring energies
on the printed and publisher HTML versions. The extraction preserves
these values without correction. The [2015 paper](https://doi.org/10.1103/PhysRevA.91.012702)
has an [NTRS catalogue record](https://ntrs.nasa.gov/citations/20160006715)
with no downloadable article. Section 3 of the 2022 article reports
poor agreement with more detailed calculations. An independent
comparison is required before using these ratios in a model. The
[Komninos and Nicolaides paper](https://hdl.handle.net/10442/6691)
identifies theoretical 2s and 2p partial cross sections near threshold
and resonance structure. It is an independent comparison lead with
separate source rights and narrower energy coverage.

Sutton's estimated errors are 30%, 20% and 10% for 4p/4s, 4d/4s and
4f/4s. His discussion on printed p. 32 warns the 4p extraction could
be biased by a factor of two or three. These relative laboratory
observables do not determine the total fraction of excited ions,
the 10 eV Maxwellian-weighted yield, or the ion recoil distribution.
The workbook and Bhatia can support I22 source cross-checks. Sutton
provides a separate I23 observable at its reported electron energy.

## Archive and literature leads with remaining gaps

- [Science Data Bank 10.57760/sciencedb.30583](https://www.scidb.cn/en/detail?dataSetId=5b6ebb2212824b28b8552af3449173e6)
  lists angle-resolved He⁺ + He capture and elastic workbooks.
  The provider's record specified `EMBARGO` at retrieval. The
  [DataCite metadata](https://doi.org/10.57760/sciencedb.30583) and
  [article](https://www.cpsjournals.cn/en/article/doi/10.1088/1674-1056/ae64d6)
  specify angles in a laboratory convention. The native normalization
  and angular axes require the actual workbook before comparison with
  the r6 He⁺ + He source. The requested I24 proton + Na channel still
  needs absolute angular bins, tails, uncertainties and producer
  provenance. Nearby IAEA entries D119989, D119991 and D119992 are
  different reactions and do not corroborate those quantities.
- [Suárez et al. (2018)](https://doi.org/10.1039/C8CP03725F) reports
  water-cation dissociation potential energy surfaces for a prepared
  excited H₂O⁺ wave packet. The RSC supplementary ZIP request
  returned HTTP 403. Its preparation and dissociation scope does
  not establish the joint outgoing photoelectron and fragment law
  over the D04 photon spectrum. The described three-body
  O⁺ + H + H branch differs from the D04 O⁺ + H₂ channel.
- [Woodruff and Samson](https://doi.org/10.1103/PhysRevA.25.848) and
  [Lindle et al.](https://doi.org/10.1103/PhysRevA.31.714) identify
  He⁺ n=2 photon partials or ratios over their reported ranges. The
  candidate Berkeley PDF responses were empty HTTP 202 responses.
  There is no new numeric import from them. Incident OH sticking
  and accommodation on lunar material remain unresolved.

## Replay and verification

```bash
python3 -B utility/chemical_database/native_m3_partial.py \
  --acquisition test_output/chemical_database/m3-source-r8-tU34BqBU \
  --output-dir /tmp/chemistry-m3-r8-new-packet
```

The recipe uses Poppler's `pdftotext` against the native text layer
and the existing lexical XLSX parser. No figures are digitized.
The recipe, helper hashes, retrieval receipts and immutable originals
identify the transformation. An altered source, missing table token
or broken importer is an error. Source reproduction does not satisfy
a production, conservation or observational gate. Formal I21–I24
remain SKIPPED pending their existing prerequisites.
