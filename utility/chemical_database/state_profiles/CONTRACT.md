# M3 source-bound physical-state profiles — implementation contract

Written before implementation, 2026-11-18, for approved P01–P09. This is
supporting infrastructure for item1; it cannot close I21–I24 or manufacture a
physical source mapping. The source/data roadmaps retain SHA-256 identities
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source copy adds approved development/resource exceptions; unchanged
scientific requirements, native metadata and qualification order control
meaning. Neither roadmap is reconciled or edited by this work.

## Inventory and entry points

Production `ReducedKinematics::ResolveState` already constructs additive
reference masses and explicit formation/excitation energies. Its current
callers construct states programmatically; there is no source-bound file
loader. `PhysicalStateProfile::Registry` will authenticate explicit profile
CSVs and source files, then call that SAME resolver. It will not calculate
new chemistry or choose a ground state, threshold, isotope mixture, PDF,
physical-error budget or speed cutoff. The independently developed physical
event composer accepts resolved states and remains separate from this loader.

The real preparation path can presently reproduce six neutral isotope masses,
four constants and separated-neutral-atom/free-electron reference-zero states
from the frozen NIST isotope/CODATA sources plus the approved P02 energy-zero
definition. No ion formation energy, molecular formation energy or excited
level is fabricated. Unsupported states remain explicit unavailable entries.

The existing mixed `reference_mass_inventory.json` also contains held ASD
thresholds. Preparation will read the two native permitted source responses
directly, not copy the mixed inventory or its held energy fields into the
published chemical database. Original research and managed-import permission
are distinct. This work never alters the immutable installed database or
canonical raw data. All derived profiles here are supporting research outputs.

## Frozen schema and refusal semantics

Native schema: `srcMoon-M3-state-profile-v1`. Caller supplies source descriptors
`id,path,sha256,role`, and exact path/SHA descriptors for three CSVs. Roles
are `isotope-mass`, `constants`, `energy`, or `reference-definition`.
All bytes must authenticate before publication of a registry. File integrity
does not establish source authority or permission; the caller separately owns
review, allowed use and formal qualification.

CSV headers, in order:

```text
isotope,neutral_mass_kg,source_id,source_record
name,value,unit,source_id,source_record
state_id,state_label,state_kind,nuclei,charge,formation_J,formation_source_id,formation_source_record,excitation_J,excitation_source_id,excitation_source_record,physical_rest_mass_kg,physical_mass_source_id,physical_mass_source_record,status,missing_reason
```

`nuclei` encodes explicit isotope/count pairs such as `1H:2;16O:1`, preserving
occurrences through integer counts. Isotope inference and natural abundances
are forbidden. Charges are integer elementary-charge units. All masses are
kg; formation/excitation energies are J. Constants retain their actual units.
State kinds are `ground`, `excited`, `unspecified`, `free-electron`.
Status is `available` or `unavailable`; unavailable entries carry a reason and
have no fabricated numerical fields. A state may be resolved only when all
required energies and provenance are present. Physical rest mass is optional,
separately source-bound, and never replaced by the additive reference mass.

The only `reference-definition` numerical use is exact zero: formation of a
separated neutral isotope atom (any explicitly stated electronic level uses
the ground atom as formation zero) or free electron, and excitation of an
explicit ground state/free electron. Molecular formation and excited-level
energy require an actual energy source. The definition references the frozen
approved decision register; it is a coordinate zero, not measured data.

Duplicate identities, malformed integers/CSV, invalid units/signs, absent
source IDs, source/hash changes and unsupported schemas are ERROR/refusal.
Missing declared source or unavailable state is SKIPPED/unavailable at the
owning runner. Failed loads/lookups preserve previous registry/output.
No missing quantity is zero, no caller tolerance is widened, and available
means only a source-bound declared reference state, not a qualified reaction.

## Test specification

| Field | Contract |
|---|---|
| ID/name | Existing M3 supporting scope; no new permanent U/I ID |
| Purpose | Actual native source/profile parsing and reference-state resolution |
| Production capability | `PhysicalStateProfile::Registry::{Load,GetState,GetPhysicalRestMass,GetConstant}` and SAME `ReducedKinematics::ResolveState` |
| Inputs | Frozen NIST native isotope/CODATA responses; approved decision register; exact profile/source SHA descriptors |
| Independent oracle | Native source strings, independent Decimal mass conversion, immutable provenance/hash/identity invariants |
| Units/frame/time | kg/J and source constant units; no trajectory frame, epoch or distribution selected by a state catalogue |
| Random controls | None; no RNG |
| Procedure | Prepare real partial profile; compile native loader probe; compare all six neutral masses/electron/constants; exercise refusal controls |
| Metrics/acceptance | Exact source identity/field semantics; relative native-load error at most twice the reported long-double machine epsilon against independent Decimal conversion; unavailable ion/molecular states stay unavailable; failed output remains unchanged |
| Expected PASS | Evaluated parser/arithmetic/provenance criteria pass |
| Failure modes | Changed source bytes, duplicate/state/units/energy-reference errors, missing state, invented physical-mass fallback |
| Required data/hashes | Actual native source and recipe hashes recorded in prepared manifest and test results |
| Predecessors | Source integrity before local kernel checks; formal M1/M2 and same-campaign gates remain unchanged |
| Artifacts | result/manifest, exact commands, native stdout, profile/source hashes, detailed comparison CSV and source audit |
| Status | PASS criterion evaluated/met; FAIL evaluated numerical mismatch; ERROR corruption/build/infrastructure; SKIPPED absent authoritative inputs |

No new physical uncertainty or activation threshold is introduced. Complete
M3 still requires source-supported channel/state/energy/PDF mappings, actual
providers/exposure/boundaries, linked tests and the frozen protocol.
