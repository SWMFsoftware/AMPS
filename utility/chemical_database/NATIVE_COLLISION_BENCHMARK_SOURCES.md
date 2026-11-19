# Original proton–helium differential benchmark source

## Source-specific admission contract — 2026-11-19

This component archives the original NASA-CR-114270 / NTRS 19700032941
report and its matching official citation metadata. It produces **no numerical
cross-section CSV**, OCR-based table, digitized curve, repaired cell or activated
reaction kernel. Original bytes and the native public-use determination control
source access; a successful download alone never grants rights or physical
qualification.

Primary source: D. H. Crandall (1970), *Differential scattering with charge
transfer to the 2s state of hydrogen from H+/He and H+/Ar collisions*,
[NTRS record 19700032941](https://ntrs.nasa.gov/citations/19700032941).
The PDF is the native scanned report, 12,816,867 bytes, SHA-256
`d0644902915fb04e9b97d8e7fca0f7342eac41e67f40ca4699909cb8d390466c`.
Its official citation JSON is 2,406 bytes, SHA-256
`d13f994bdb99ff31e200914097c0b641c52b62007c24fce95ff7dafca1c3fb6e`.
Exact original HTTPS request/final URIs, HTTP 200 and actual retrieval UTC are
retained in the new packet and acquisition evidence.

## Rights are attached to this record only

The narrow reviewed policy is `ntrs-19700032941-public-use`, host
`ntrs.nasa.gov`, restricted to paths identifying this exact citation. The
actual citation declares `distribution=PUBLIC`,
`disseminated=DOCUMENT_AND_METADATA`, and
`copyright.determinationType=GOV_PUBLIC_USE_PERMITTED`; its copyright-marker
and third-party-material indicators are false. The scanned first page states:
“This document has been approved for public release and sale.” This notice
was visually inspected in the actual PDF, rather than accepted from OCR.
Both the citation and original notice remain available for review.

This is a source-specific public-use determination, **not an assertion that
this contractor thesis is a US Government-authored public-domain work**, a
Creative Commons licence, or an application of the unrelated NASA CEA Apache
licence. All native ownership flags remain visible in the complete record;
false legacy flags are not silently reinterpreted as independent proof of
ownership. The [official STI API dictionary](https://sti.nasa.gov/docs/OpenAPI-Data-Dictionary-062021.pdf)
identifies the metadata fields; that reference link is not imported as another
cleared numerical payload. General NASA guidance explicitly distinguishes
public release from unrestricted reuse of copyrighted contractor works;
this operation relies on the record's explicit public-use determination and
native notice, not NASA hosting alone. No external contact occurred.

## Measured quantity scope and unresolved interpretation

The report contains H+–He measurements over 3–20 keV and selected laboratory
angles roughly 1°–7°. This overall study range is not a complete rectangular
energy-angle support. Individual tables and real blanks control available
measurements; no outside-domain extrapolation or all-angle normalization follows.

| Native table / printed page | Source quantity / units | Native error / interpretation condition |
|---|---|---|
| 8 / 165 | Total differential scattering, 10^-16 cm²/sr; columns indexed by theta*T in keV degrees | Estimated absolute error ±15%, except near 1° ±30%; a printed decimal ambiguity remains unresolved. This is total scattering, not pure charge exchange. |
| 9 / 166 | Source charge-transfer probability P0, dimensionless; retain its native definition | Estimated absolute error ±10%; no new confidence level, covariance or independent-state partition is inferred. |
| 10 / 167 | H(2s)-selected charge-transfer probability; native column/scientific-notation scaling remains in the PDF | Estimated absolute error ±40%; not a ground-H or all-final-state total. |
| 11 / 168 | Differential H(2s) capture, 10^-18 cm²/sr; theta*T-indexed source columns | Estimated absolute error magnitude ±40–50%, not an asymmetric +40/-50 interval or assigned Gaussian sigma; native blank cells remain missing. |
| 12 / 169 | Explicit 6.25 keV angle case; separate total-scattering, P0, P2s and H(2s)-capture quantities | Each quantity retains its own source scale/error. No one table-wide uncertainty or shared independent-error assumption is imposed. |

These are source descriptors, not typed numerical cross-section tables.
The native report covers additional processes, including Ar, whose physical
interpretation has not been promoted by this narrow He-table review.
NTRS machine-readable OCR is not imported as a numerical table. The ambiguous
printed cell is not repaired from an expected curve or guessed decimal.

## Reproduction, API and status

The new evidence directory is:

```text
test_output/chemical_database/physical-model-integration-20261119T075920Z/collisions/
```

Its `source_packet.json` contains the original PDF and native citation JSON,
with direct source identities and pinned rights evidence. No shared parser,
API, source inventory or canonical file is changed by this note. The shared
catalogue's source-access API can expose the archived original by source ID;
there is no numerical `table-access` descriptor or physics evaluator here.
Installed applications resolve archived objects without following the original
research paths.

Source integrity checks compare expected size and SHA before copying into a
new output, verify copied bytes, and retain original acquisition receipts.
The maintained packet loader additionally enforces the reviewed policy and
record/path boundary. Missing source means SKIPPED; changed bytes or malformed
metadata mean ERROR. A source-admission PASS means the declared original is
available under this record policy; it is never a U/I benchmark comparison,
complete CEX daughter model, source-state mapping or M3 milestone PASS.
