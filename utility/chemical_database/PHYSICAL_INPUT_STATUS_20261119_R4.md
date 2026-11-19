# Additional physical inputs — r4 source disposition

This is the database-first, primary-publication follow-up requested on
2026-11-19. Source HEAD is `26f1a65c58a8a70f1c6301b1a6c4d283a4558c94` plus
the recorded working-source hashes. Its separate immutable release is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r4`.
Select that directory's own `database-lock.json`; the original root lock and
earlier releases are not retargeted. Source access is not physical activation,
linked verification, observational validation or M3 closure.

The pre-change scope, entry points and acceptance contract are in
`PHYSICAL_INPUT_EXTENSION_CONTRACT.md`. Native originals, retrieval receipts,
packet selection and detailed research reports are retained under
`test_output/chemical_database/m3-physical-inputs-20261119T231417Z/` in the
source checkout. The installed release retains native admitted payloads,
metadata, exact selected packets, declared recipes and helper dependencies;
querying it does not require that historical checkout path.

## New numerical inputs and what they establish

| New source | Native data supplied | Scientific limits preserved |
|---|---|---|
| [ExoMol OH+ MoLLIST](https://www.exomol.com/data/molecules/OH_p/16O-1H_p/MoLLIST/) | 815 actual state rows and 12,044 Einstein-A transitions for the source-labelled isotopologue, with exact definition and native units | Definition declares 823 states. All actual transition references resolve, but missing rows are not created and completeness is not inferred. Relative spectroscopic energies do not supply the absolute neutral–ion formation offset. |
| [ExoMol OH MoLLIST](https://www.exomol.com/data/molecules/OH/16O-1H/MoLLIST-OH/) | 1,878 actual states and 54,276 transitions, independently sourced from the earlier MYTHOS product | Definition declares 2,640 states; optional native quantum columns disagree with their JSON labels. Only the four mandatory state fields are semantically named; later tokens remain positional. No guessed v/e/f reassignment. |
| [ExoMol H2 ARLR](https://www.exomol.com/data/molecules/H2/1H2/ARLR/) | 2,051 states and 32,915 transitions | Dataset-relative internal energies and radiative coefficients are source inputs, not a PHIDRATES formation/Q/state-population crosswalk. No omitted-exit bound or unprovided coefficient uncertainty is invented. |
| [ExoMol H2 RACPPK](https://www.exomol.com/data/molecules/H2/1H2/RACPPK/) | 302 actual states and 4,712 transitions | Definition declares 303 states. This alternative H2 calculation is not merged with ARLR or normalized as an extra channel. |
| [IAEA CollisionDB proton–Na search](https://db-amdis.org/collisiondb/search/?reactant1=H%2B&reactant2=Na&process_types=52) | Seven previously absent excited-Na AOCC records D107579/580/582–586; 49 exact printed coefficient rows and native pages | These are capture-plus-ionization composites, not pure CEX or H(n)-resolved partials. The equation is not rendered and the upper energy endpoint is absent. Native 10% model experience error is not experimental 1σ uncertainty. No unsupported fit evaluator is activated. |
| [Zenodo 5044539](https://zenodo.org/records/5044539) | Original Apollo 15221 hydrogen/water TPD and energy/site-distribution author files, four exact lexical tables | Keep measured TPD separate from model components and inverted site/energy distributions; retain signed values and preparation/specimen assumptions. Neither these nor earlier Apollo samples define a universal lunar sticking law. |
| [Zenodo 14588634](https://zenodo.org/records/14588634) | Six author XLSX originals, eleven worksheets and 19,151 represented cells for soil water adsorption/migration | Native cell coordinates, strings, gaps and metadata survive. Source-file versus final-paper figure-number drift is unresolved. Isotherm/migration measurements are not silently converted to prompt sticking, accommodation or a time-dependent lunar reservoir law. |
| [Zhao et al., Nature Chemistry 2021](https://doi.org/10.1038/s41557-021-00730-1) | Three publisher source-data ZIPs for OH(A)+H2 quenching, including selected reproduced experimental distributions and calculated rotational/vibrational/wavepacket products; thirteen lossless long-cell views | Preserve measured/model distinctions, source missing markers, duplicate headers, blank separators, collision energies and Ntot=0. These fixed-energy/state data are not a thermal total rate, complete joint product law or an electronic-quenching rule for all environments. Do not renormalize the two displayed reactive outcomes while omitting elastic/inelastic survival. |

ExoMol sources use their pinned CC BY-SA4.0 notice; IAEA records use the pinned
CollisionDB CC BY4.0 notice; both Zenodo deposits explicitly declare CC BY4.0.
The OH-quenching supplementary grant is tied to its own primary article's
CC BY4.0 evidence. These grants do not clear other providers' independent
datasets. Exact packet/table/source IDs, units and original hashes are queried
through the same installed Python, CLI and C++ Catalogue API as earlier inputs.

## Correction to the earlier MYTHOS interpretation

The newly acquired native MYTHOS `.def` version20240526 explicitly labels
`131898.260953 cm^-1` as **Higher energy with complete set of transitions**
(line64). Its SHA-256 is
`4da22ac7d223b6b92174f716e52c80ecf4fa77ca8387d21f56c457f455816f4c`.
The earlier JSON field called the quantity `max_energy`. Comparing it with
the largest retained state energy, `132196.65537 cm^-1`, as a metadata defect
was incorrect: a transition-completeness energy is not a maximum state energy.

`amend_exomol_metadata.py` reproduces the r4-only correction from the exact
old packet and new version-matched native evidence. The amended packet keeps
all twelve original source IDs and all native/CSV/recipe member bytes, removes
the misleading equality flag and records the native label and provenance.
Earlier release bytes and historical QA remain unchanged. The maintained
`native_exomol.py` no longer infers a state cutoff from the ambiguous JSON key.
This correction is distinct from the real MoLLIST/RACPPK declared-count gaps.

## Publication and database searches that did not supply admitted data

The atomic follow-up retains metadata for the [CODATA2018 author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC9888147/)
and the [CODATA2022 primary paper](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf).
The latter was published in2025 (JPCRD54,033105,
[DOI10.1063/5.0279860](https://doi.org/10.1063/5.0279860)); it is not a2024
paper. Its explicit rights-reserved notice and missing source-specific managed
reuse grant are preserved. Generic NIST/PMC access is not treated as permission
for these numerical payloads. No indexed snippet or manually copied threshold
is used as database science data. A failed EuropePMC XML request remains an
acquisition ERROR, not a source record.

ATcT requests returned access challenges; challenge bodies are not
thermochemistry. Additional WebBook and primary-record pages retain their
actual source/permission state. Three native LAMDA OH collisional files provide
rotational relaxation information, with managed-operation permission unresolved;
they remain metadata-only. Rotational OH(X) data are not OH(A) electronic
quenching, and O ground-term fine structure is not O(1D)/O(1S) relaxation.
OPEN-ADAS still requires written permission for managed insertion. Existing
ASD, CCC, LXCat and other source-specific holds are not relaxed.

Following IAEA records into [Hill et al.2023](https://doi.org/10.1088/1741-4326/acf5da)
identifies Na collision references including
[Abdurakhmanov et al.2021](https://doi.org/10.1103/PhysRevA.104.042820).
The publication lineage is useful research evidence, but it is not a newly
acquired complete differential array or a proven original citation for
D119990. Do not conflate Na+H ion-pair thermal rates with proton+Na capture.

## Remaining admission and qualification work

New radiative matrices reduce missing-source gaps. A physical molecular
provider still needs an isotope/state/energy-zero/formation-threshold crosswalk,
source-compatible physical masses, declared complete exits or omitted-branch
bounds, initial populations and source-supported uncertainty/model-error scope.
The source's relative energy zero alone cannot define these quantities.

The approved heavy-ion reservoir approximation does not require a complete
joint outgoing-electron PDF; it does require source-consistent incoming states,
rate/domain/moment evidence, signed charge/momentum/energy ledgers and the
approved numerical/model-error checks. None of the new composite Na records
supplies separated CEX products or a forward-angle error bound. Existing P07
Na3s reference approval remains distinct from recovered experimental provenance.

The approved ideal absorbing collector remains a legitimate controlled
verification boundary. Real lunar ion retention, neutralization, accommodation,
readsorption and thermal release require their separate supported material,
temperature, incidence, charge/state and distribution domains. A TPD inversion,
reflected-charge conditional law or sputter yield cannot substitute for all
incident-particle outcomes.

Full linked convergence/statistical matrices, global budgets and formal M1/M2
predecessors remain required. No reaction is enabled, uncertainty assigned,
observation selected or acceptance threshold changed by this release. Exact
source/software test counts and new release identities are in its work-package
results and publication/installation receipts. A source/API PASS establishes
that stated operation only; M3 remains EXPERIMENTAL / NOT QUALIFIED.
