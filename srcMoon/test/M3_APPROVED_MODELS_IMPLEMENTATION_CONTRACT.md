# Approved P01–P09 implementation contract

2026-11-18, before implementation. Source revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus existing local development.
The user approved all nine recommendations in
[the decision proposal](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md).
Approval covers the explicit reduced-model, collector and laboratory scopes;
it supplies no absent state/threshold/joint-kinetic/material data and does
not waive M1/M2, U24, D06 or any formal I21–I24 requirement.

## Authority and inventory

Source roadmap SHA256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
External roadmap SHA256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds approved development-order/resource amendments; these govern
work order. Common formal gates, native metadata and data-root rules still
control qualification/interpretation. Neither roadmap is changed here.

`Photochemistry` supplies the complete D04 scalar families and source means,
with physical daughter support false. `ReactionExecution` supplies number,
integer material, reference-mass, momentum and energy plan checks;
`MoonReactionAdapter` commits actual PIC packets. `MoonReactionController`
admits constant histories and explicit providers only. Both force tracers
and guarded common chronology exist. Neither a physical event model nor
an ion-material kernel can be inferred from these mechanical services.

The existing exact PIC/catalog mass/charge checks are retained. New reference
mass configurations must use actual consistent compiled tables; the reviewed
legacy Na configuration/reference remains unchanged. Source-only math is
linked into the normal Moon library, with thin independent fixtures testing
those same kernels. No new shared-core physics or generated-copy patch is
planned. Canonical validation data remain read-only.

## Implementation and source dependencies

1. Record P01–P09 as approved, retain all activation requirements, and acquire
   source-bound isotope/constants/state evidence where publicly available.
   A reproducible extraction must distinguish isotope masses, formation
   energy, state/threshold data and source uncertainty. Missing molecular
   thresholds/states remain missing. No source-dependent state is invented.
2. Implement P02/P03/P04 in dependency-free production `ReducedKinematics`:
   additive reference masses, explicit energy zeros, reviewed error limits,
   parent-rest photon/two-body recoil and mass-weighted microcanonical breakup.
   Inputs carry source identity and explicit states; no internal RNG, implicit
   ground state, energy truncation or unsupported full-family activation.
3. Implement P08 explicit strict/ideal-collector runtime choice and a
   source/state-resolver installation guard. Actual sphere/terrain dispatch
   must precede neutral sticking and preserve actual impact geometry/time.
   The artificial reservoir records incident material/charge/momentum/energy;
   it implies no neutralization, heat, implantation depth, sputtering/release
   or surface-charging feedback. Missing state/energy information is refused.
4. Implement P09's conditional laboratory He/silica law only within the
   publication's actual regimes/velocity definition. Reproducible TableI
   processing and source QA precede any laboratory capability claim. An
   unavailable primary convention means an explicit unavailable branch.
5. Implement the P07 Na partial-channel native-node rate and experimental
   forward-COM products with physical Q, retaining the existing total-rate
   registries. Source metadata, not a total-fit shortcut, determines scope.
6. P05 explicit state/cascade support and P06 signed reservoir interfaces
   retain their missing-data guards. Any physical provider that cannot supply
   every requested state/energy/bath term remains disabled. The first kernel
   may use its declared parent-rest frame; this is not a completed moving
   LSO photon provider. Document exactly which runtime adapters exist.

## Supporting test specification

| Field | Contract |
|---|---|
| ID / name | Supporting extensions within U15/U17 and existing chemistry/boundary meanings; no new permanent U/I IDs. |
| Purpose | Verify approved mathematical closures, exact accounting, explicit configuration and production dispatch without declaring physical validation. |
| Production capability | Source-bound mass/state helpers, reduced recoil, bounded partial CEX, strict/collector boundary and conditional He laboratory operator. |
| Entry points | New species kernels linked through normal Moon makefile; existing conservative planner/PIC adapter and actual sphere/terrain callbacks where operational. |
| Configuration | Explicit opt-in models, source identities, states, frames, energy-zero/validity bounds; unchanged baseline by default. |
| Inputs | Declared analytical fixtures plus exact native source tables with hashes; supplied open uniforms/normal variates, physical event/impact state and explicit bath/collector. |
| Independent oracle | Closed-form two-body identities; independent Jacobi-coordinate distribution; exact integer identities; independent source fit/node evaluation; actual owned gas/collector inventory. |
| Units | SI kg,m,s,m/s,C,J; source eV,u,cm²,keV/u conversions explicitly sourced. |
| Frame/time | Parent-rest inertial recoil scope distinguished from LSO; actual SPICE ET and moving facet rest frame for production impacts. |
| Deterministic/stochastic controls | Explicit draws and fixed independent original-parent seeds; actual1/8MPI where linked; no inferred independence between daughter records. |
| Procedure | Contract/source QA, standalone production-kernel checks, registered supporting children, isolated clean build/current M0/nightly and matched native dispatch controls. |
| Metrics | Nuclei/charge/reference-mass/energy/momentum closure; source/domain/error refusal; two-body partition; Jacobi moments/distribution; actual callback collector transfers and baseline preservation. |
| Acceptance | Existing U/I criteria unchanged. Representation-level mathematical bounds are frozen before checks; derived sampler distribution is the oracle, never current output. No new physical error allocation is selected implicitly. |
| Expected PASS | Only the actually evaluated declared software/laboratory criterion passes. Full IDs remain SKIPPED where physical data or predecessors are absent. |
| Failure modes | Wrong mass scaling, lost H multiplicity/electron, implicit state, mean energy confused with PDF, partial/total CEX substitution, source-gap interpolation, neutral sticking used for ions, false lunar retention claim. |
| Required data | Frozen D04/D05 for existing scope; native NIST/constants/state and IAEA/Tellez sources for new source-bound cases; no observations generated. |
| Data hashes | Recorded before use in source manifests and each runner; no fabricated hash for unavailable source data. |
| Predecessors | Actual same-campaign formal requirements unchanged. Early supporting development remains explicitly permitted. |
| Artifacts | stdout/stderr, result/manifest/input/source/build/data identities, actual commands, independent CSV/PNG/EPS comparisons where applicable, budget/collector records and resource histories. |
| Status semantics | PASS evaluated/met; FAIL evaluated criterion miss; ERROR build/parser/corruption/runtime/infrastructure; SKIPPED absent required data/capability/configuration/predecessor. |

## Execution and preservation

Use the existing global supporting protocol, child exits0/1/2/77 and retained
runner/log/artifact paths. Source-only tests do not configure AMPS. Native
packages use an isolated checkout outside the source tree, remove `build`
before normal Config.pl/SPICE/make, and share the approved3600-second wall
and200GB aggregate-RAM limits, including setup/plots/restoration/cleanup.
Default8MPI/OMP1, plus explicit1-rank controls. Preserve the original failed
attempts and reviewed four-rank nightly reference/zero-byte diff. Longer
statistical campaigns are deferred explicitly rather than reducing gates.

No active-root Makefile/site change, canonical data write, staging, commit,
push, external message or automatic reference freeze is authorized. STATUS
and component/application/test READMEs must report actual completed versus
unavailable runtime paths and exact retained evidence.
