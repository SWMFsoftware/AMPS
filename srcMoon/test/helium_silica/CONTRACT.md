# P09 He+/silica conditional charge-survival contract

Written before implementation, 2026-11-18. P09 is a user-approved development
choice, not a new U/I ID or an observational-validation result. The source
roadmap SHA712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17
and external roadmap SHAd10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4
differ by approved development-order/status amendments. They retain their
physical/qualification gates. Canonical data/raw observations remain read-only.

| Field | Declared contract |
|---|---|
| ID/name | Supporting P09 source-specific laboratory kernel; no new permanent U/I ID. |
| Purpose | Reproduce Tellez2012 Eq3 and all seven TableI source fits without turning a conditional backscattered charge fraction into an incident implantation/retention probability. |
| Production capability | New dependency-free species/HeliumSilicaChargeSurvival registry and source-coordinate evaluator. Normal library integration belongs to parent; no boundary dispatcher is enabled here. |
| Entry points | FrozenV1Sources, Registry::LoadPublishedTable/GetFit/SelectFit/EvaluatePublishedCoordinate; physical-velocity conversion explicitly unavailable pending source convention evidence. |
| Configuration | Exact versioned source PDF and deterministic processed TableI CSV identities; source normal-incidence145degree silica experiment only; no arbitrary material/geometry extrapolation. |
| Inputs | Target Si,16O,18O or source O(mean); incident energy[eV]; positive mechanism-specific reciprocal-velocity sum[s/m] supplied in the exact source coordinate. No energy-to-velocity/elastic-scattering model is inferred. |
| Oracle | Actual TableI extracted reproducibly from native PDF text, visually independently checked against its rendered table; independent high-precision Decimal evaluation of the published exponential. |
| Units | Vc[m/s]; inverse velocity[s/m]; incident energy[eV]; dimensionless logP/P. TableI Vc values/uncertainties multiplied by exactly1e5. |
| Frame/time | Laboratory surface-rest frame,4He+ incident normal to isotope-exchanged silica;145degree scattered detector. Source-coordinate evaluation has no inferred epoch or lunar material/rotation coupling. |
| Stochastic controls | None; deterministic source law. No random charge assignment, BCA event or experimental observation is generated. |
| Procedure | Verify PDF identity; native-text extract/QA seven rows; pin resulting CSV; compile actual production kernel; compare complete table metadata and independent exponential across all rows; check endpoint gaps, source corruption, invalid inputs and transactional refusals. |
| Metrics | Source/table identity; exact row/parameter/domain mapping; probability/logP error; preserved no-fit gaps; output immutability; reported uncertainty endpoints. |
| Acceptance | Every source row/number/domain matches the reviewed native table; numerical law/envelope within64binary64eps relative to independent high-precision arithmetic for positive representable outputs. No physical accuracy or observational agreement tolerance is introduced. |
| Expected PASS | Evaluated local formula/transcription controls pass; physical velocity mapping and independent laboratory agreement stay unavailable/unqualified. |
| Failure modes | Missing source, changed hashes, wrong1e5 conversion, using fit intercept as a probability prefactor, extrapolation through0.8–2keV O gap, mixing normal/full speed conventions, probability->retention conversion, assigning a confidence level to source± values. |
| Required data | Actual acquired Tellez etal APL101151602 DOI10.1063/1.4758699 publisher PDF, TableI and Eq3; no digitization/OCR. |
| Data hashes | PDFbf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2; deterministic processed table identity recorded after extraction before kernel tests. |
| Predecessors | This supports source-specific arithmetic only. Normal Moon boundary/ion fate, formal M1/M2 and chemistry/surface predecessor gates remain unchanged. |
| Artifacts | Source processing command/tool/version/hash manifest, actual tableCSV/QA, compile/stdout/stderr, source/compiler/executable identities, resultJSON, independent reference/candidateCSV+PNG+EPS. |
| Status | PASS evaluated local criterion; FAIL evaluated mismatch; ERROR corrupt source/build/parse/receipt/contradiction; SKIPPED missing source/dependency. Unavailable physical mapping is not a failed or passed laboratory validation. |

## Scope and domain decisions

Source TableI gives CIN Si Ei<=7.0keV, O CIN Ei>2.0keV, and O AN
Ei<0.8keV. Intersect these with the actual study600–7000eV. Thus O800eV
and2000eV are unavailable, along with the entire intervening range. The
caption mentions the600–800eV AN estimate, but the row inequality remains
strict; no endpoint extension is inferred. Study endpoints600/7000eV are
retained where the row inequality includes them. Si has no published AN row.
O(mean) is a separately published coefficient, not an isotope-abundance law.

Eq3 is P+=exp[-Vc*(1/vi+1/vf)]. Intercept a belongs to corrected-yield Eq4
and is NOT a prefactor in P+. Reported Vc± maps monotonically to an envelope;
the source does not define a confidence level/covariance here. It is not
presented as a1sigma or95percent interval. a/r²/O:Si values are preserved as
metadata; no independent uncertainty or observation is invented.

The paper says reciprocal-velocity definitions depend on mechanism and cites
Brongersma2007. General primary theory distinguishes nonlocal perpendicular
AN and local atom-relative velocities. That does not, without explicit source
mapping, establish which coordinate the numerical TableI fit used. Therefore
only the explicit source-coordinate law is available. A raw-vector/energy
adapter must refuse while its exact fitted-coordinate mapping is unresolved;
no normal/full-speed or elastic outgoing-energy assumption is silently made.

This is conditional P(He remains/returns+1 | the selected backscattered class)
including survival/reionization, not P(retained | incidentHe+). It supplies no
reflection probability, isotope mixture, outgoing angles, implantation depth,
surface electron/heat ledger, or lunar residence/re-emission. Source fits alone
and formula verification do not validate their predictions against independent
measured spectra. Full lunar ion-surface behavior remains disabled.
