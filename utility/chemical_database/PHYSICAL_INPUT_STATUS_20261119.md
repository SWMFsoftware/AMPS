# Physical-input status — 2026-11-19

**Later same-day research update:** this document describes the r2 acquisition
and remains its historical source disposition. The subsequent repository report
`srcMoon/test/M3_PHYSICS_LITERATURE_20261119.md` supplies a source-derived Na+
mass, actual CHIANTI O I radiative arrays, new molecular/Apollo measurements
and concrete reduced-model options. Some earlier missing-input statements now
need source/model admission rather than a new measurement. The installed r2
snapshot and its original documentation remain immutable; no provider is
activated by the later research.

**Approved r3 integration:** the subsequent immutable
`snapshots/amps-chemistry-20261119-r3` release imports the reviewed CHIANTI O I,
Apollo TPD, Chang water-product, JWST and ARTEMIS author files, with reversible
or exact lexical views. Application-independent model definitions and explicit
mathematical/reduced C++ helpers are supplied separately. Missing metadata,
physical-model admission, managed-source holds and application qualification
remain explicit; neither r1/r2 nor canonical sources are edited. See the new
installation README section9, source-specific notes and
`PHYSICAL_MODELS_INTEGRATION_CONTRACT.md`. Historical findings below retain
their date/scope and are not reinterpreted as current universal data absence.

This is the maintained disposition of acquired source inputs and remaining
physical requirements for the lunar M3 work. It records the source tree at
`8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129` plus the working changes.
It activates no physical process, changes no immutable release or canonical
validation file, and closes no milestone. Source packets authenticate bytes,
provenance, permission and native quantities. **Source access and lexical table
access are separate from a physically qualified provider and linked U/I PASS.**
The present acquisition package adds material; it is not evidence that every
remaining needed datum has been found or is globally unavailable.

## 1. Acquired inputs and their remaining scope

| Input / primary source | Actually acquired and usable source scope | Remaining physical or permission condition |
|---|---|---|
| [NASA CEA](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da) | Pinned `thermo.inp`, native compiler-format authority, Apache 2.0 LICENSE and data-path NOTICE; full 2,111-record source, with exact projections of 19 selected gas species and their coefficient arrays. Native temperature/enthalpy conventions and original row citations remain attached. | Thermochemical model coefficients are not isotope-resolved physical rest masses or resolved electronic/rovibrational populations. Assigned H(0) is not automatically formation Hf(0); electron reference and thermal increments must be retained. Uncertainties/covariance and a PHIDRATES threshold crosswalk are not supplied. The complete source has eleven nonincreasing condensed-phase intervals, preserved as failures; the selected gas intervals are increasing. No evaluator may silently repair or use an invalid interval. |
| [ExoAtom Na / Na+](https://www.exomol.com/exoatom/) | Native NIST-derived ExoAtom 20240901 definitions and spectrum-relative state lists: **423 Na I states and 186 Na II states**, energy and source uncertainty in cm^-1. ExoMol's separately published product is acquired under its [CC BY-SA 4.0 policy](https://www.exomol.com/data/licence/), retaining the underlying NIST lineage. | Both native mass fields equal **22.98976928 Da**, and both `Ionisation` fields are **null**. Separate E=0 spectrum references do not supply the neutral–ion formation offset, physical ion rest mass, isotope shifts, lifetimes or D04 threshold compatibility. Never derive a physical Na+ mass by accepting the equal source masses. This product-specific grant does not change the hold on independently acquired NIST ASD files. |
| [ExoMol H2O DTU](https://www.exomol.com/data/data-types/photo/H2O/1H2-16O/DTU/) | Newly acquired separate **423.15 K and 573.15 K** absorption files, **7,933 rows each**, native support 110.35–230.364 nm, and exact lexical projections. Source isotopologue, temperature and cm²/molecule units remain attached; managed storage allowed under CC BY-SA 4.0. | Native files contain **679 and 1,448 negative samples**, respectively. Preserve signed source values: they are not an admitted nonnegative probability/cross-section provider, and no clipping, smoothing or repair is authorized. Absorption is not channel yield, solar rate or daughter kinetics. No colder-temperature or missing EUV extrapolation. |
| [ExoMol OH MYTHOS photo product](https://www.exomol.com/data/data-types/photo/OH/16O-1H/photo-MYTHOS/) | Previously retained 300 K total-cross-section file plus photo definition version **20241206**, cm²/molecule and nm; byte access allowed under CC BY-SA 4.0. | Actual **19,173 rows** differ from metadata **19,172**. This older photo-product discrepancy is separate from the new rovibronic acquisition below. Do not drop an endpoint. Absorption does not give photofragment state, recoil or cascade branching. |
| [ExoMol OH MYTHOS rovibronic product](https://www.exomol.com/data/molecules/OH/16O-1H/MYTHOS/) | Newly acquired native definition version **20240526**, compressed state list with **10,606 states**, and **10,000-row partition table**, plus exact lexical projections. Native E and uncertainty in cm^-1, lifetime field in s and dimensionless partition Q remain source quantities. | The actual state maximum **132196.65537 cm^-1** differs from declared **131898.260953 cm^-1**. Preserve both; no cutoff or corrected definition is inferred. States/partition data alone do not supply molecular formation energies, complete radiative exits or a full daughter law. A supplied lifetime is not a complete state-to-state branching matrix. |
| [IAEA CollisionDB](https://db-amdis.org/collisiondb/) | New 21 state/shell-labelled e–Na/e–He integral-fit pages and exact lexical projections: **172 native printed coefficient rows**. IAEA's [data licence](https://db-amdis.org/collisiondb/licence/) grants CC BY 4.0 absent a dataset-specific exception; original attribution is retained. Existing D119990 partial-capture knots and distinct D107581 line/patch sources remain separate. | Native exact exports failed HTTP 500, so printed precision is not promoted to full-precision source data. Some rendered formula/domain information is absent or internally inconsistent; no repaired evaluator or tail is implied. Alternative overlapping fits must not be summed. No relevant joint electron/ion daughter distribution was acquired. A patch is plot geometry; a composite cross section is not pure CEX. |
| Frozen [D04 PHIDRATES](https://phidrates.space.swri.edu/) and [D05 Voronov](https://www.pa.uky.edu/~verner/col.html) | Existing immutable canonical rates remain the approved rate baseline: four atomic, seven H2O and four OH_exp reactive channels, distinct quiet/active 1 AU endpoints and native mean energies; D05 supplies the selected scalar electron-impact coefficients. | Existing canonical use is preserved. A new managed numerical copy is not cleared by the current PHIDRATES policy; Voronov research-evaluation permission is distinct from shared managed-storage permission. Means and scalar Maxwellian rates do not establish joint daughter laws. D04 source thresholds must not be replaced by coarse fit thresholds from another source. |
| [TU Wien fixed Apollo 16 sputter deposit](https://researchdata.tuwien.ac.at/records/hjn8d-hyh84) | Explicit all-data CC BY 4.0 grant; native archive size/MD5/SHA and exact 18-table reproduction pass. **Four measured tables / 17 rows**, separate from **14 model-output tables / 2,033 rows**. Incidence angle, amu/ion, atoms/ion and native errors remain unchanged. | Target sputter yield is not incident reflected/implanted/retained fraction or outgoing charge distribution. Atoms/ion uses the authors' bulk-stoichiometry conversion. The [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12267057/) says sample 68501 in abstract/results and 68051 in methods; retain the unresolved crosswalk. No source temperature dependence or error confidence level is inferred. |
| [Tellez 2012 He+/silica](https://doi.org/10.1063/1.4758699) | Existing source-specific published-coordinate conditional charge-survival law and Table I processing support a bounded laboratory benchmark. | Conditional charge among selected backscattered particles is not unconditional neutralization or lunar retention. The exact full-versus-normal reciprocal-velocity coordinate tied to Table I is unresolved. Shared managed-table rights remain unknown; research evaluation and numerical redistribution are different operations. |
| Existing [NIST O I ASD](https://www.nist.gov/pml/atomic-spectra-database) packet | Native levels/lines and selected upper-state records are retained with source uncertainties and missing A values. | Direct ASD managed-storage permission remains unresolved. Complete radiative exits, missing coefficients, isotope/state mapping, quenching/polarization and any omitted branch bounds must be established before a total cascade lifetime or unit-normalized decay law. ExoAtom access does not clear these independent ASD records. |

The new H2O temperature files and OH rovibronic/partition records are in
`test_output/chemical_database/physical-input-integration-20261119T062501Z/exomol-native/source_packet_v2.json`.
They are distinct from the older photo products and Leiden inspections in
`test_output/srcMoon/m3-molecular-databases-20261118T083551Z/`. The older
HDF5 negative-sample/partial-closure findings remain in force; its initial
report's lack of HDF5 tooling does not override subsequent actual inspection.

## 2. Numerical acquisitions held from managed storage

| Source | Actual source evidence and hold | What could change the disposition |
|---|---|---|
| [CCC atomic sources](https://atom.curtin.edu.au/CCC-WWW/) | A 900-member He archive (899 integral tables plus README), Na `TICS.3S`, and a proton–He source file were acquired. Native README calls the He results unpublished as of 2026-11-09. Eight nonincreasing energy pairs in one table remain unsorted. No explicit managed-data grant was located: metadata-only. | Actual source-specific licence/permission covering managed numerical storage, followed by state/domain/unit QA. These integral totals still do not become joint product distributions. |
| [OPEN-ADAS](https://open.adas.ac.uk/terms-and-conditions) | Current terms retain the explicit written-permission requirement for managed-database insertion. ADF01/ADF07 numerical payloads stay held; no contact or permission was obtained. | Explicit written permission for the intended operation. Do not relabel a private cache as an exception. |
| [SpuBase v1.0](https://zenodo.org/records/10783295), [v2.0](https://zenodo.org/records/18741496), [v2.1](https://zenodo.org/records/22139055) | All three native archives pass producer MD5. Zenodo metadata says CC BY 4.0, while the embedded source licence says GPLv3. This scope conflict is preserved: metadata-only. Each version has 72 mineral fit tables; 3D fit columns in v2.1 occur only in 28 SW files. | An unambiguous version/file-level grant or reviewed resolution of mixed licence scope. Do not choose the more permissive notice silently or treat simulation curves as measurements. |
| [Leiden](https://home.strw.leidenuniv.nl/~ewine/photo/), [MPI-Mainz](https://www.uv-vis-spectral-atlas-mainz.org/uvvis/), [UMIST](https://umistdatabase.uk/downloads), [KIDA](https://kida.astrochem-tools.org/help) | Useful retained spectra, aggregate partials and network metadata, with source-specific issues. Current managed-storage permission remains unknown. KIDA modelling-download evidence is retained separately from a shared managed-database grant. | A source-specific grant for the new managed operation; then a declared source/field/domain contract. ISRF network rates must not be used as 1 AU solar rates. Leiden's published solar normalization is explicit, not a universal factor for another network. |

Permission decisions are operation-specific, not declarations that every local
research use is forbidden. See repository policy `utility/chemical_database/RIGHTS.md`. A public URL, successful
HTTP response, open article or another provider's licence does not clear a
numerical dataset. No permission request or author contact was performed.

## 3. Remaining physical provider requirements

1. **State and energy mapping:** bind each actual parent/product isotope,
   electronic/rovibrational state, physical rest mass, additive reference mass,
   formation/excitation energy, source threshold and uncertainty to compatible
   source definitions. The prior Na D04-versus-ASD difference of about
   **5.98895654e-6 eV** is not resolved by the new ExoAtom relative lists or
   equal mass fields. Neither rounding provenance nor an uncertainty interval
   is invented to admit it. Keep source and numerical errors distinct.
2. **Photon products:** all 15 D04 reactive hazards and distinct source endpoints
   remain represented. Mean excess energy alone supplies neither a full
   conditional energy distribution nor angular/polarization/state correlations.
   Molecular channels need actual state/threshold mappings; an unqualified O
   label is not automatically O(3P). Absorption totals are not extra deaths.
3. **Electron and charge-transfer products:** a prescribed Maxwellian rate
   is not an electron-sharing/recoil distribution. State-dependent beam totals
   require an explicit relative-energy convention, supported domain and bath
   population; missing high-energy tails cannot be set to zero silently.
   Full CEX additionally needs complete final-state branches, physical Q and
   differential energy-angle support or an independently bounded reduced law.
4. **Ion surface:** first impact, ideal absorption, selected backscatter charge
   survival, sputtering and long-term retention are distinct observables.
   A lunar law needs material/temperature/incident state and energy-angle
   domains plus normalized reflected/neutralized/implanted outcomes and their
   emitted distributions. No new acquired table supplies all those quantities.
5. **Cascade and molecular relaxation:** missing radiative branches, OH/H2
   internal-state relaxation and environment-specific quenching cannot be
   omitted without a source-supported bound. A selected subset of A values
   cannot be normalized to a complete lifetime by assumption.
6. **Linked qualification:** correct event assembly still requires real
   trajectory exposure, common chronology, source/bath accounting, exact
   committed-particle checks, convergence/statistics and frozen numerical
   criteria. I23's measured forcing dependency and I24/U24 applicable surface
   behavior are not waived by local tests or a source-access release.

## 4. Approved bounded closures that remain possible

The user-approved P01–P09 decisions in repository document
`srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md` authorize explicitly named reduced experiments; they do not supply missing
numbers or make approximations physically validated.

| Approved bounded scope | Required actual source/error evidence before use | Honest claim if its checks pass |
|---|---|---|
| Complete D04 scalar rate/selection controls | Authenticated original archive, endpoint/units/state labels, all hazards and existing source-domain refusal | Local rate/selection verification; no daughter kinetics or M3 closure |
| Fixed source-mean isotropic two-body photon closure | A resolved source/state/Q/mass crosswalk, physical photon/bath ledger and actual model-error bounds | Conservative named reduced event model; no measured joint fragment law |
| Mass-weighted microcanonical three-neutral-body closure | Explicit O/cofragment states and thresholds, source-supported available energy and approved error assessment | Statistical reduced H2O breakup; not dissociative-ionization kinematics |
| Na ground-product partial CEX at native measured nodes | The actual partial-channel identity, native supported incident energies, states/Q, source uncertainties and explicit forward-CM approximation bound | Partial-channel reference transport; not total CEX or unqualified solar-wind pickup |
| Ideal absorbing collector | Actual incident material/charge/momentum/energy ledgers and a declared mathematical boundary | Verification and first-impact diagnostic; not implantation, neutralization, sticking or release |
| He+/silica published-coordinate benchmark | Actual Table I/source identity, fit domain and uncertainties; leave vector projection unavailable | Conditional laboratory charge-survival operator; not lunar retention |
| TU Wien sputter-table reference comparison | Retained specimen mismatch, exact measured/model roles, native energy/incidence and error conventions | Source-specific sputter-yield comparison; not an incident-ion fate kernel |

The last row describes a newly available source-specific comparison opportunity,
not automatic approval of a production sputtering model. New model choices,
unsupported sources and any changed physical interpretation still require
review. Model and numerical error bounds must be actual and independently
checked; the approved statistical/numerical allocation is not a licence to
widen a failed physical source gate.

## 5. Immutable source identities and reproducibility

The active manifests below are selected explicitly by a new release build,
never edits to a published immutable release. Paths are historical acquisition
and reproduction locators in the AMPS repository, not installed runtime
dependencies. The portable catalogue resolves its own archived objects; it
does not require this checkout or follow original paths to recover data. NASA
`atomic/source_packet_native.json` is an earlier native-only version of the
same dataset ID: select the current `atomic/source_packet.json`, not both.

| Packet relative to the AMPS checkout | Source records | SHA-256 at this report |
|---|---:|---|
| `test_output/srcMoon/chemical-input-acquisition-20261119/atomic/source_packet.json` | 7 | `57e309a23da8bda56cfea214b6d31d3685d4fc8f90f83051a5f9fe2566f62d29` |
| `test_output/srcMoon/chemical-input-acquisition-20261119/atomic/source_packet_exoatom.json` | 6 | `23599c0b7c572f2435f212dee7f3ae013392daa69bafaa5bfe72936da9ba3e8b` |
| `test_output/srcMoon/chemical-input-acquisition-20261119/collisions/source_packet.json` | 48 | `e8b4220e9fef205027797be7a606cc2aecb8e438223eecbcb9455a4d0928eb31` |
| `test_output/srcMoon/chemical-input-acquisition-20261119/surface/acquisition-20261119T062501Z/source_packet.json` | 20 | `e3120e6b8206eba08062b74743d9aba7c7a1ee52bc6c90aa59819ae0537f9c5a` |
| `test_output/srcMoon/chemical-input-acquisition-20261119/surface/acquisition-20261119T062501Z/spubase_source_packet.json` | 6 | `3e421264d47f1a93173a4ad9f0d7794f2677bcb693e85418a720ae1f1f36eedf` |
| `test_output/chemical_database/physical-input-integration-20261119T062501Z/exomol-native/source_packet_v2.json` | 12 | `f777959c00f8bfda6c5297dcb5e0255469b423ef1ca25e994756d78a9380139d` |

Detailed acquisition receipts and scientific caveats are in the adjacent
research packets, especially the surface `REPORT.md`, collision `REPORT.md`,
`atomic/processed/nasa_thermo_native_index.json`, and the molecular packet's
`REPORT.md` plus subsequent HDF5 inspections. Each source packet pins exact
HTTPS retrievals, UTC, hashes, native parent identity, rights evidence and the
maintained extraction recipe where a derived table exists. No numerical
value, ordering, gap, negative sample or uncertainty is repaired for import.

The new TU Wien exact-member source recipe passed nine real-source controls.
New IAEA lexical projections preserve all 172 printed rows. These results
establish only the declared source/software criteria. Missing source metadata,
unsupported physical capability and permission-held numerical data remain
SKIPPED / NOT VALIDATED; corruption and parser/runtime failure remain ERROR.
An evaluated physical mismatch remains FAIL.

Unsuccessful ATcT challenge responses, IAEA export HTTP 500 responses, NILS
posterior HTTP 404, instrument-directory HTTP 401 and publisher HTTP 403 are
preserved research outcomes. They establish what failed in these specific
bounded attempts, **not proof that no authoritative data exist elsewhere**.
No challenge/login body is a physical source record.

Formal execution still follows the existing roadmap predecessors. Approved
early implementation does not turn a missing M1/M2 predecessor into PASS or
permit observational scoring before the frozen operator/data/selection gates.
