# Source permissions and the inventory contract

This file records the source-specific policy used by `source_inventory.py`.
The user authorized implementation and a local installation at
`/data/vtenishe/CHEMICAL_DATABASE`. That authorization is not a third-party
licence. This inventory neither copies scientific files nor enables physics.
It reports every known dataset, including an explicit missing-packet entry.

## Contract written before implementation

`inventory(amps_root: Path, validation_root: Path) -> list[dict]` performs
read-only enumeration. Dataset records contain stable IDs, source/version
descriptions, independent storage/evaluation/redistribution permissions,
evidence, capabilities, issues and native file records. File records retain
absolute original paths, native logical names, SHA-256 identities, media
types and roles. Optional table records describe original column names,
native units and actual record counts; they do not rewrite numerical values.

The permission enum is `allowed | pending | prohibited | unknown |
not-applicable`. **Only `allowed` authorizes numerical managed storage.**
`metadata-only` is the resulting availability, not a permission enum.
Unknown sources remain discoverable, with no numerical payload copied.
The builder must enforce the file-level permission as well as its dataset
policy; an allowed provider does not clear mixed-source reports/manifests.

The inventory is deterministic for unchanged inputs, uses no network or RNG,
performs no normalization or fitting, and writes no source/data file. A new
research packet is not READY merely because it parses or its licence allows
use. `source-access` and `table-access` mean lexical access only, not a
rate law, complete reaction family, event sampler or physical qualification.
The independent software checks are actual native CSV/ASCII counts/headers,
unique stable IDs, strict ADAS holds, permitted-source table coverage and
explicit missing-packet records. They are not observational tests.

## Operation-specific decisions, checked 2026-11-18 and 2026-11-19

| Source | Local managed numerical storage | Basis and restrictions |
|---|---|---|
| IAEA CollisionDB | allowed | [Official licence](https://db-amdis.org/collisiondb/licence/) explicitly applies CC BY 4.0 to data and metadata unless a dataset states otherwise. Preserve original authors, dataset URL/DOI when supplied, licence and modifications. The retained D119990 page states no dataset-specific exception; its generic website copyright footer does not replace the data licence. This clears storage, not missing initial-state/bibliographic/uncertainty or full-kinetic metadata. |
| NASA CEA thermochemical data | allowed for the pinned NASA deposit | [NASA CEA](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da) includes Apache 2.0 LICENSE.txt and NOTICE.txt explicitly identifying data/thermo.inp and data/trans.inp. Retain both notices and the exact commit. Only the acquired thermo.inp is imported here. This does not grant a different NIST compilation's rights or turn thermochemical reference masses/enthalpies into an isotope/state rest-mass profile. |
| CHIANTI native O I | allowed with attribution | The pinned [official use notice](https://www.chiantidatabase.org/) and [database referencing guidance](https://www.chiantidatabase.org/referencing.html) cover use and incorporation with version/original-reference attribution. Preserve exact native files and references. No rights are transferred to direct NIST ASD exports or mixed derived mass records. |
| Chang2021 supplementary water-product workbooks | allowed for matching CC-BY4.0 articles | Each supplementary ZIP is tied to its own article DOI and pinned publisher/EuropePMC XML with the explicit article CC-BY4.0 statement. DOI10.1038/s41467-021-22824-7 does not clear DOI10.1038/s41467-021-26599-9; both grants are retained separately. Keep source labels, fits and data distinct. |
| Jones2020 Apollo TPD, JWST OH and Poppe2022 pickup deposits | allowed for the individually pinned CC-BY4.0 records | Native Zenodo records3822365,10634184,6994648 supply the file-specific producer identities and licence. Preserve original author/source attribution. The grant for a deposit does not clear a neighbouring publisher paper or turn modeled columns into measurements. |
| AMPS authored model definitions/helpers | AMPS software under Apache2.0 | Original mathematical software and factual source descriptors retain the repository LICENSE.txt. This software licence does not clear cited third-party numerical data, select parameters or activate physical models. No acquisition record is fabricated for locally authored files. |
| NASA-CR114270 / NTRS19700032941 | allowed for this original PDF and official record | The exact [official citation record](https://ntrs.nasa.gov/citations/19700032941) declares public distribution and `GOV_PUBLIC_USE_PERMITTED`, with no third-party-material/ownership marker; the original cover approves public release and sale. This narrow record-specific determination is retained with the PDF. It is neither a blanket public-domain claim for NASA contractors nor a licence for OCR-derived science tables. |
| TU Wien lunar-soil sputtering deposit | allowed for version hjn8d-hyh84 | [Native record](https://researchdata.tuwien.ac.at/records/hjn8d-hyh84), DOI10.48436/hjn8d-hyh84, explicitly applies CC BY 4.0 to all its data files. Import the pinned archive/member bytes with attribution, distinguishing measured and simulated tables. This does not provide incoming-ion trapping or charge-state distributions. |
| ExoAtom Na/Na+ source files | allowed as the separately published ExoMol product | [ExoMol data licence](https://www.exomol.com/data/licence/) covers all ExoMol data under CC BY-SA 4.0. Preserve the ExoAtom product's explicit NIST lineage and version; do not apply this grant to direct ASD downloads. The native definitions' missing ionization field and equal neutral/ion mass remain limitations. |
| SpuBase Zenodo versions | unknown: metadata-only | Native record licences indicate CC BY while acquired embedded notices assert GPLv3. Both remain evidence; the conflict is not silently resolved by choosing the more convenient grant. The numerical ZIP payloads are held. |
| CCC He 2026 native collision archive | unknown: metadata-only | The acquired native README marks the calculations unpublished (November 2026), without an established managed-storage grant. Availability is not permission. No payload is admitted under another collision provider's licence. |
| ExoMol | allowed | [Official licence](https://www.exomol.com/data/licence/) applies **CC BY-SA 4.0**, not CC BY 4.0, to all ExoMol data. Retain attribution, licence link, change record and relevant articles. Sharing adapted material has ShareAlike conditions. Local use does not resolve the OH row-count discrepancy or turn absorption into channel kinetics. |
| Author-deposited enstatite / Apollo17 EELS data | allowed for the acquired deposit objects | Native Zenodo records [15978325](https://zenodo.org/records/15978325) and [10794539](https://zenodo.org/records/10794539) explicitly give `metadata.license.id=cc-by-4.0`. Keep attribution and distinct simulation/experimental roles. The full EELS archive was not acquired; this does not authorize copying neighbouring publisher PDFs. |
| NIST unmarked constants/isotope source exports | allowed, under the stated public-information policy | [NIST Copyrights & Disclaimers](https://www.nist.gov/copyrights-disclaimers), “Use of NIST information,” permits copying/distribution except for material marked copyrighted. The retained CODATA ASCII and isotope response and their respective [constants](https://pml.nist.gov/cuu/Constants/index.html) / [isotope](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses) source pages have no product copyright marker. Preserve NIST/source citations and the actual policy; do not label all NIST data public domain. |
| NIST ASD / ground-state threshold exports | unknown: metadata-only | The [ASD homepage](https://www.nist.gov/pml/atomic-spectra-database) explicitly marks the product copyrighted, with all rights reserved and possible future charges. [NIST's SRD policy](https://www.nist.gov/open/license) distinguishes SRD from its non-SRD grant. No product-specific managed-database permission was established in this audit. An unmarked individual CSV/CGI response does not erase the product notice. O I radiative sources and mixed derived inventories containing ASD thresholds therefore stay on hold. |
| Other NIST atomic/electron research | unknown: metadata-only | Do not transfer the unmarked constants/isotope decision to SRD107 or another NIST compilation. Product-specific permission was not established here. |
| OPEN-ADAS | pending: metadata-only | Native ADF01/ADF07 terms explicitly require written permission for insertion into a managed database and redistribution with modelling code. See [official terms](https://open.adas.ac.uk/terms-and-conditions) and retained `adas_he_cex_2.dat` lines86–96. A private path/cache is not an exception. No provider was contacted and no permission was obtained. |
| LXCat | unknown: metadata-only | The retained terms restrict third-party redistribution/commercial use; [redistribution policy](https://nl.lxcat.net/instructions/redistribution.php). Source/contributor-specific local managed-storage permission was not established. The challenged numerical export did not succeed. |
| PHIDRATES | unknown: metadata-only | The [source](https://phidrates.space.swri.edu/) describes its scientific calculation service and carries a copyright notice. No explicit managed-storage licence was located. Existing canonical D04 bytes and legacy use remain untouched; this inventory does not silently authorize a new managed copy. |
| Voronov / UKY | unknown for managed storage; research evaluation allowed by the source statement | The [official introduction](https://www.pa.uky.edu/~verner/info.html) explicitly discusses research/code use and requests acknowledgement; [source table/program](https://www.pa.uky.edu/~verner/col.html). This is not a located blanket managed-database/redistribution licence. Preserve that distinction rather than applying another provider's policy. |
| Leiden photo tables | unknown for managed storage | The [official site](https://home.strw.leidenuniv.nl/~ewine/photo/) supplies numerical research data and citation instructions but no explicit managed-storage grant was located. Native negative samples and incompatible versions remain genuine issues. |
| MPI-Mainz | unknown: metadata-only | [Atlas](https://www.uv-vis-spectral-atlas-mainz.org/uvvis/), [linked imprint](https://www.mpic.de/3626961/imprint). No source-specific managed-data grant was established. An open-access article licence must not silently be applied to every externally linked data file. |
| UMIST RATE22 | unknown: metadata-only | [Official downloads](https://umistdatabase.uk/downloads) provides the network; an explicit dataset grant was not located. The captured New Mills copyright is webpage software, not the RATE22 data licence. |
| KIDA | unknown for managed storage | [Help](https://kida.astrochem-tools.org/help) explicitly permits downloading reaction sets/storing reaction lists for modelling, while [legal notices](https://kida.astrochem-tools.org/legal) reserve site intellectual property. Preserve this useful research-use evidence; do not silently upgrade it to a cleared shared managed-data licence. |
| Tellez silica table / publisher PDF | unknown for managed numerical storage: metadata-only | The retained publisher-PDF cover permits private study/research copying but restricts further distribution; [original](https://pure.tue.nl/ws/portalfiles/portal/3503988/583864115491353.pdf), DOI10.1063/1.4758699. Paper, extracted table and future release rights are separate. Local P09 mathematical verification is not a licence grant for a shared managed database. |
| Other papers, Phelps/Rice and material deposits | unknown unless individually evidenced | Keep bibliography and availability; do not infer a numerical database licence from a paper's accessibility, a hosting site's unrelated licence, or another dataset's grant. |

“Unknown” records an evidence gap; it is not a legal claim that every local
research use is prohibited. The approved strict inventory policy holds the
new managed operation until its scope is established. Installation should
not touch or disable existing canonical sources or legacy application modes.

## File grouping and scientific caveats

Provider-specific files within research packets are split by their original
logical names. Mixed reports, inspectors and manifests can embed numbers
from several providers, including ADAS; those files stay `unknown` rather
than inheriting IAEA/ExoMol permission. The same rule holds for the derived
NIST mass inventory containing both unmarked constants/isotopes and ASD
ionization thresholds. Raw constants/isotopes are separate records.

Expected byte hashes follow a declared provenance order: complete packet or
raw/acquisition inventories, then source-processing provenance, then the
current native inventory, then retained exploratory metadata. An obsolete
exploratory descriptor cannot retarget a current output by reusing its old
basename. Its differing expectation remains a historical-metadata issue;
no hash is selected because it happens to agree with the observed bytes.

Native IAEA D119990 CSV nodes are an allowed lexical table with the exact
`energy_eV_per_u,cross_section_cm2` header and twelve rows. D107581 line and
patch arrays remain distinct; the latter is plotting geometry and neither
is a pure-capture total. ExoMol `.photo` tables retain their headerless native
ASCII bytes; documented source columns/units are metadata, not a new header.
No value is converted, clipped, gap-filled, sorted or renormalized here.

Permission is independent of scientific interpretation. Missing A values,
source states, covariance, domain bounds, complete radiative exits or joint
daughter laws remain missing. Existing source-integrity PASS results are
not promoted to model verification or validation. No licence request,
author contact, production build or scientific comparison is performed by
this inventory.

## Targeted remaining-input additions, 2026-11-20

`figshare-13274840-cc-by-4.0` covers only Laffon's sticking Table1 deposit
13274840.v1/file25563521. The exact native API record declares CC BY4.0 and
lists the payload download URL. Admission checks the article ID, version,
licence URL and actual file membership. This does not clear another Figshare
deposit. Retain the source DOI/attribution, gold-relative XPS normalization,
20K substrate/300K gas conditions, native uncertainties and substrate distinctions.
The undeclared native eight-bit text encoding is not silently converted.

`nist-jres-121005-public-use` is specific to the inspected Burgess H/O
thermochemical evaluation, DOI10.6028/jres.121.005, its official PDF/PMC7339710
rendition and the pinned NIST technical-publication notice. That notice
explicitly includes the Journal of Research and specifies reuse/attribution
for NIST-employee work. Retain the requested NIST acknowledgement and article
citation with derived views. This is not a grant for the separately copyrighted
WebBook/ASD compilations, another NIST paper or a source published elsewhere.
Table cells preserve their native state, zero-K, reference and uncertainty
semantics; permission is not a physical source-to-reaction mapping.

## November 24 archive inputs

- `zenodo-4555716-cc-by-4.0` is tied to the pinned native deposit record,
  its listed `diff_cross_sections.zip` and the producer size and MD5.
  Keep the Pawel Swaczyna attribution, data DOI, CC BY 4.0 notice and
  change record for the two CSV views. The approval covers the deposit and
  its CSV views. The journal paper and independent datasets have their
  own licences. [Zenodo record](https://doi.org/10.5281/zenodo.4555716)
- The live D119990 HTML is admitted with the existing IAEA policy and
  exact numeric comparison to the retained sodium partial. The HTML
  and numerical projection have separate identities.
  [IAEA policy](https://db-amdis.org/collisiondb/licence/)
- The NORAD He photoionization text files are catalogued with their
  native identities and the OSU CC BY-NC 4.0 notice. The general shared
  snapshot retains metadata for these files. OSU identifies its contact
  for commercial permissions. [NORAD licence notice](https://norad.astronomy.osu.edu/)
- Curtin CCC, NIFS, Verner–Yakovlev, NIST and the cited experimental
  publications remain metadata-only. Preserve the terms of the source
  objects individually. The NILS public counts carry CC BY-NC-SA 4.0
  and the Bodi paper carries CC BY-NC 3.0.
  [NILS record](https://doi.org/10.5878/v4sy-3t95)
  [Bodi publisher notice](https://doi.org/10.1039/C4SC00791C)

## Source-specific publication follow-up

The `nist-jres-water-1966-kim-1992` policy applies to the exact reviewed
Dibeler and Kim technical articles and their lexical CSV views. Preserve
the original citation and the acknowledgment “Republished courtesy of
the National Institute of Standards and Technology.” SRD 107 and the
separate original experiments cited by Kim require their own review.
[NIST technical-publication notice](https://www.nist.gov/open/license)

The `mdpi-atoms13090074-cc-by-4.0` policy covers Cui et al.'s exact
published article and the Table 1 projection. The grant binds the
publisher PDF and its embedded CC BY 4.0 notice. Cite the DOI, authors
and the lexical conversion. Angular data requested separately from the
authors require their own file-specific review. [Cui et al.](https://doi.org/10.3390/atoms13090074)

The r8 policies bind managed bytes to the exact reviewed versions and
records listed in [the source report](NATIVE_M3_PARTIAL_SOURCES.md).
`figshare-12389750-v2-cc-by-4.0` applies to version 2, file 22860461,
with CC BY 4.0 attribution to Kalha et al. and Yeh and Lindau.
`ntrs-19720018903-public-use` applies to the exact Sutton report and
records NASA's source-specific public-use determination.
`mdpi-atoms10010009-cc-by-4.0` applies to the inspected Bhatia article
and the projections of his own calculations. The Science Data Bank He
dataset and Suárez publication record enter this snapshot as metadata.

The r10 `iop-aceff9-cc-by-4.0` policy binds Wang and Bruggeman's exact
IOP PDF, its CC BY 4.0 statement and the corresponding Crossref version
of record. Credit the title, authors, journal and DOI. The new
`ntrs-19730006011-public-use` grant is specific to Sutton and Kay's
NASA preprint and the matching public-use determination. The three
CSV projections identify their source columns and preserve the
authors' uncertainties. [Source-specific rights](NATIVE_M3_LAB_BENCHMARKS_SOURCES.md)

The Moore–Mapaye Table 1 and Laurell Extended Data Fig. 4 permissions
are tied to the exact files, DOI and CC BY 4.0 notice for each article.
The Crowley evaluation is CC BY 3.0. Its underlying primary experiment
has separate rights. The new IAEA member is covered by CollisionDB's
CC BY 4.0 terms. The NIFS report, Shingal calculations and helium
line-emission articles retain catalogue descriptions and individually
attributed facts. Their source files have separate access policies.
[Photon](NATIVE_M3_PHOTON_BENCHMARK_SOURCES.md)
[Ice](NATIVE_M3_OH_ICE_SOURCES.md) [Sodium](NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md)
[Helium](NATIVE_M3_HELIUM_EMISSION_SOURCES.md)

The `aas-peschel-2026-cc-by-4.0` policy covers the identified Peschel
et al. article and the attributed projection of its Table 1 model inputs.
The article carries its own CC BY 4.0 notice on page 1. Its citations have
their own access and reuse terms. [Article](https://doi.org/10.3847/PSJ/ae4901)

## I22 decision-2 sources, 2026-11-26

`rsc-bodi-2014-cc-by-nc-3.0` covers only Bodi et al. 2014, Chem. Sci. 5,
3057, doi:10.1039/C4SC00791C. It applies to the RSC article PDF as served by
the PSI DORA repository, SHA-256 `8882498b…`, and to the projection of its
printed OH⁺ appearance energies. Every page of the article prints "This
article is licensed under a Creative Commons Attribution-NonCommercial 3.0
Unported Licence", so the article is its own notice.

The redistribution field is `non-commercial-with-attribution`. Use is
non-commercial only, and any redistribution must be non-commercial with
attribution. This is the first non-commercial grant in the catalogue. It does
not extend to the article's supplementary files or to the works it cites.

The other records in this packet are metadata-only:
- **Schörghofer 2023 and 2025.** Their CC BY 4.0 payloads await a genuine
  HTTP-200 retrieval; iopscience returned 403.
- **OSTI, HAL and GUP repository copies.** No licence statement was located.
- **Publisher PDFs supplied by the user.** They are identified by SHA-256 in
  record metadata only.

The user's approval to incorporate the sources follows this policy. It is not
treated as a third-party licence. [Sources](NATIVE_I22_DECISION2_SOURCES.md)
