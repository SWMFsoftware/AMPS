# P05 O I radiative-source acquisition and QA

This source contract was written before implementing the processing script on
2026-11-18, after the user approved P01–P09. It assigns no new permanent U/I
ID and activates no physical decay model. The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the external roadmap is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source roadmap includes approved early development and resource/status
amendments. Native NIST metadata controls interpretation; unchanged formal
qualification and predecessor requirements continue to control validation.

| Field | Contract |
|---|---|
| ID/name | Supporting P05 NIST O I levels and radiative-transition source acquisition. |
| Purpose | Obtain the actual public O I exports without an arbitrary wavelength window, preserve their provenance and gaps, and audit listed exits from the 2s2.2p4 1D J=2 and 1S J=0 states. |
| Production capability/entry | Source preparation only. No production kernel, PIC state, parser, decay law, lifetime, collector or chemistry dispatcher is modified or activated. |
| Configuration | Exact official NIST ASD request parameters, raw response bytes and source version retained. Include allowed and forbidden transitions and missing transition probabilities; apply no minimum strength or accuracy threshold. |
| Inputs | Official NIST level and line responses, query forms and interpretation help. All downloaded files belong in a new immutable-within-this-package research directory under test_output. |
| Independent reference | Native levels/line identities, source-unit metadata, source line counts, and cross-check of selected exits against a separate unrestricted query for the exact upper-level identity where supported. Source bibliographic references are preserved. |
| Units | Retain source energy/wavelength units and Einstein A units exactly; parse numeric values only where source metadata establishes the units and notation. Preserve accuracy codes instead of inventing statistical uncertainties. |
| Frame/time | Atomic rest-frame spectroscopy, not a lunar observation epoch or surface state. Retrieval UTC and database version are provenance, not a source event time. |
| Stochastic controls | None. No generated observations, fit, digitization or parameter tuning. |
| Procedure | Preserve complete queries and native responses; verify successful response format and species; parse native records reproducibly; bind levels to upper/lower line states using source identities; preserve blanks/qualifiers; enumerate every listed exit of each selected state; audit completeness separately from a successful download/parse. |
| Metrics/acceptance | Raw/processed SHA-256 identities, native and parsed row counts, unique level mappings, selected transition counts, missing/ambiguous A values, duplicate component identities, bibliographic references, and any explicit source completeness limitation. A source-integrity PASS means only the evaluated integrity checks passed. |
| Expected PASS | Reproducible source capture/parse/identity checks pass. A complete physical lifetime is not declared merely because every returned A value is numeric. |
| Failure modes | HTTP error/HTML gate, empty or truncated export, changed schema, wrong ion/species, ambiguous upper level, missing A, unresolved component combination, arbitrary wavelength/strength cut, selection after model residuals, or unsupported completeness inference. |
| Required data/hashes | Official NIST ASD O I response bytes. Actual request/response and output hashes are recorded at retrieval/processing; no invented source identities are prefilled. |
| Predecessors | No model qualification is inferred. Physical cascade activation still requires reviewed complete applicable transition data, source-state mapping, quenching/alignment/fine-structure treatment and the existing M1/M2/M3 predecessor gates. |
| Artifacts | Raw native exports/forms/help, request/header/retrieval manifests, source-processing script identity, normalized CSV records, selected-exit CSV, QA/provenance JSON and a source-gap report. |
| Status | PASS for a criterion actually evaluated; FAIL for an evaluated mismatch; ERROR for request/parser/corruption/infrastructure failure; SKIPPED/UNAVAILABLE for missing authoritative data or completeness metadata. No missing probability is set to zero or a unit branch. |

The source package is bounded to one hour including retrieval, processing, QA
and reporting. Canonical validation data and raw observations remain read-only.
No total lifetime or normalized branching probabilities are synthesized while
the applicable radiative exits are incomplete or not authoritatively bounded.

Before adding descendant processing, the unrestricted native export was
inspected for the fine-structure states reached by the selected upper levels.
The source audit also retains the stored graph reachable from those levels,
including distinct 3P J states and every listed multipole component. Graph
closure means closure of the retrieved stored records only. A node without a
stored exit is not assigned an infinite lifetime or zero physical decay rate.
