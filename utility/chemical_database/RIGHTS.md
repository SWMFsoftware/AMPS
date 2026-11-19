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

## Operation-specific decisions, checked 2026-11-18

| Source | Local managed numerical storage | Basis and restrictions |
|---|---|---|
| IAEA CollisionDB | allowed | [Official licence](https://db-amdis.org/collisiondb/licence/) explicitly applies CC BY 4.0 to data and metadata unless a dataset states otherwise. Preserve original authors, dataset URL/DOI when supplied, licence and modifications. The retained D119990 page states no dataset-specific exception; its generic website copyright footer does not replace the data licence. This clears storage, not missing initial-state/bibliographic/uncertainty or full-kinetic metadata. |
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
