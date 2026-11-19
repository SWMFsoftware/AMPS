# CHIANTI atomic source tables and explicit reduced-model scope

The optional CHIANTI source packet makes native atomic levels and radiative
quantities available through the application-independent chemical catalogue.
It does not select a Moon species/provider or qualify a radiative/chemical model.
Use the explicit release/profile lock and table IDs
`chianti-oi-levels-five` and `chianti-oi-radiative-lowest-five`.

The [official O I page](https://db.chiantidatabase.org/o/o_1.html) attributes
observed energies to NIST (downloaded2010), theoretical energies to Zatsarinny
and Tayal (2003, ApJS148,575), and transition probabilities to Froese Fischer
and Tachiev (2004, ADNDT87,1, DOI10.1016/j.adt.2004.02.001). The original
`o_1.elvlc` and `o_1.wgfa` remain complete, with seven levels/sixteen entries.
The selected projection retains the lowest five explicitly indexed levels and
all nine original links, including six exits from 1D2/1S0 and three ground
fine-structure links. Weak 639.2/295.8-nm exits are preserved, not renormalized
away; the green branch produces 1D2 that can subsequently decay.

The acquired ion-description page says11.0.1; the download page advertises
database11.0.2 and a separate IDL11.0.5. These labels must not be conflated.
The packet pins actual native SolarSoft bytes and metadata; it does not claim a
release-tar identity that was not checked. Data footers retain original credits.

[CHIANTI's homepage](https://www.chiantidatabase.org/) states that the package
is freely available and asks appropriate acknowledgement. Its
[referencing guidance](https://www.chiantidatabase.org/referencing.html)
explicitly addresses incorporation into other databases: retain the version
and original calculation citations. The named `chianti-use-attribution` policy
uses these archived operation-specific statements. No SPDX identifier, grant
for unrelated ASD products, or licence for the original papers is invented.

## Native conventions and missingness

`native_chianti.parse_ion_bytes(levels,radiative)` implements the fixed-column
current CHIANTI8+ elvlc/wgfa format described by technical reports
[1](https://www.chiantidatabase.org/tech_reports/01_elvlc/chianti_report_1.pdf)
and [10](https://www.chiantidatabase.org/tech_reports/10_wgfa/chianti_report_10.pdf).
The immutable returned Catalogue owns its source strings; returned mappings
are defensive copies. Malformed/nonfinite/duplicate/unresolved data refuses
transactionally. An optional explicit `maximum_bytes` checks combined input
size and never truncates science rows. An invalid native energy ordering is
preserved as a source-quality issue, not repaired. No interpolation or fit occurs.

Observed and theoretical energies remain separate **cm^-1** columns, with the
native **-1** missing marker explicit. No theoretical fallback is applied.
Wavelengths retain their original signed **vacuum Angstrom** values; a negative
value means an unobserved/theoretical endpoint. Zero wavelength is a special
two-photon/autoionization-style rate, not a single photon with inferred energy.
The gf column is dimensionless and the decay column is **s^-1**. The files do
not provide per-entry A/energy uncertainties; printed precision is not an error
bar. A whole native file is not evidence of a complete higher-level cascade,
collisional quenching, external radiation field or an observational operator.

`levels_csv(catalogue,indices)` preserves explicit requested level order.
`transitions_csv(catalogue,indices)` preserves original transition order and
refuses selection of an upper state while omitting a known downward endpoint.
It produces lexical values and identities, without normalized branch fractions.
Other applications can reuse the parser for this supported format and supply
their own explicit index selections, source policy and physical model admission.

The current selected A-set gives radiative sums corresponding to approximately
116.22s for1D and0.74695s for1S, and green/UV photon-number ratio16.0673.
These are derived theoretical/community-model values, not measurements.
Modern simultaneous Mars observations give16.5±0.4 (weighted-mean standard
error; calibration uncertainty separate) and a larger sample gives15.8±2.4
(background-noise uncertainty), supporting a near16 ratio. The older terrestrial
nightglow9.8±1.0 remains a distinct source. Never mix an observed branch ratio,
unrelated lifetime and another source's A-values without explicitly defining a
new model. [Gérard2020](https://doi.org/10.1038/s41550-020-1123-2),
[Soret2023](https://doi.org/10.1029/2023JE007762),
[Slanger2006](https://doi.org/10.1029/2006JA011972).

## Sodium candidate: reusable derivation, source admission still explicit

For a declared neutral isotope and its singly ionized electronic ground state,
CODATA's mass relation is
`m_ion = Ar(neutral)*m_u - m_e + IE_J/c^2`. Every term requires an authenticated
source, units, energy zero and uncertainty convention. A consistent reference
bookkeeping can choose ion reference mass `m_neutral-m_e` and ion formation
energy `IE_J`; physical rest mass then includes that formation contribution.
This definition is a software/model candidate, not an assigned database state.

The optional packet deliberately contains **no derived numerical Na+ mass or
copied Na ASD threshold**: the current direct-ASD managed-source hold must not
be cleared by CHIANTI's permission for its distributed O I files. The research
candidate descriptor records required source roles/formula/uncertainty method
without importing held numerical inputs. Constants/isotope sources already
permitted elsewhere remain separately discoverable; a complete numeric profile
requires independent admission of the threshold source or a newly reviewed
alternative. Physical missingness and operation-specific permission are separate.

## Verification

```
python3 -B test/chemical_database/test_native_chianti.py --output-dir NEW
python3 -B test/chemical_database/test_native_chianti.py \
  --levels PACKET/raw/o_1.elvlc --radiative PACKET/raw/o_1.wgfa \
  --output-dir NEW_WITH_SOURCE
```

The first command evaluates labelled parser/ownership/domain fixtures and
honestly skips the real-source class. The second authenticates the actual
originals and independently checks whitespace-native tokens, source-row hashes,
all weak exits and exact Fraction sums/ratios. Logs, source identities and
PASS/FAIL/ERROR/SKIPPED counts are preserved. This is shared software/source QA,
not an M3 physical or observational qualification and not a new U/I ID.
