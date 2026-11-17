# M2 applied initial thermal state — contract before implementation

The audit found a production/qualification mismatch: a generic document could
claim periodic thermal qualification while each of the 288 numerical cases
actually initialized uniform columns. A statement about the epoch alone did
not identify the applied thermal state. This extension changes no numerical
or observational acceptance limit and creates no permanent test ID.

| Test specification field | Contract |
|---|---|
| ID/purpose | Existing I19, with I20 consuming its evidence; prove all numerical cases apply the explicitly qualified physical initial state. |
| Production entries | Normal `InitializeProductionState` and its input-selected existing periodic spin-up; read-only native initial-state exporter before source/weight evaluation and cohort seeding. |
| Configuration | `thermal_initialization_mode=declared-uniform` remains the diagnostic default; opt-in `periodic-spinup` applies the existing production solver. Periodic cycle end must equal the declared simulation start; no interval is invented. |
| Inputs | Reviewed state qualification, actual pilot input/manifest/receipt/full-state CSV, current executable/source identity and canonical physical-input identity. |
| Oracle | Equality to the reviewed **applied initial condition**, not independent physical truth. Scientific method/thermal qualification requires its separate retained evidence. |
| Units/frame/time | K; metres at explicit layer centres; logical body-fixed surface elements; native TDB seconds plus exact input UTC spelling. |
| Stochastic control | Initial thermal state is before seed-dependent particles; all 288 selected cases must retain identical qualified initial-state bytes. |
| Procedure | Validate frozen review and method evidence, actual pilot application, explicit initialization method/epoch and full CSV shape; verify each actual staged input and native state before reading numerical differences. |
| Metrics/acceptance | Exact qualified state SHA-256, complete cell/layer membership, positive finite K, consistent layer depths, method/epoch/physical identity and source/executable identity. Existing thermal convergence gates remain unchanged. |
| Failure modes | Unapplied periodic claim, missing native exporter, uniform state mislabelled periodic, case-specific state/input change, incomplete or duplicate grid, wrong depth/time/method, altered metadata/data bytes. |
| Data/hashes | Model-only initial-condition evidence with path/SHA-256 descriptors; never observations, calibration data or validation truth. Raw external data remain immutable. |
| Predecessors | All existing same-campaign I13–I18, data, source, numerical and statistical gates remain required. |
| Artifacts | Native `receipt.txt.I19_initial_thermal_state.csv`, native initialization metrics, pilot/source/method evidence, per-case state and physics identity, formal aggregate membership. |
| Status | Missing exact qualification/native capability is SKIPPED; malformed, corrupt or contradictory evidence is ERROR; existing evaluated numerical/physical criteria retain FAIL. No automatic approval or READY artifact is generated. |

## Explicit qualification document

The statistical-support document's existing
`actual_epoch_thermal_qualification` object must now have schema
`srcMoon-M2-initial-thermal-state-qualification-v1`. The earlier generic
READY/evidence-list object is unavailable evidence, not permission to infer
an applied state. Required fields are:

* `status=READY`, `review.approved=true`, nonempty `review.reviewer`,
  timezone-qualified `review.approved_utc`, and
  `review.frozen_before_numerical_campaign=true`;
* `initialization_mode` (exactly `declared-uniform` or `periodic-spinup`),
  `simulation_start_utc`, `ephemeris_time_tdb_s`,
  `physics_identity_sha256`, `executable_sha256`, and
  `production_source_signature_sha256`;
* `pilot_input`, `pilot_manifest`, `pilot_receipt`, `native_state`: explicit
  path/SHA-256 descriptors of real retained pilot artifacts;
* `state_shape.logical_surface_element_count` and `state_shape.layer_count`;
* `method_qualification.status=READY`, explicit `criteria` and nonempty
  `evidence_files` descriptors, together with separate nonempty top-level
  `source_evidence_files` descriptors. Uniform initialization is not banned:
  its qualification must explicitly apply to **that uniform state**, rather
  than attaching evidence for an unexecuted periodic state.

Descriptor paths are absolute or resolved against the statistical-support
document. The native pilot manifest binds its actual input and binary. The
qualification's physics identity is recomputed by the shared canonical
projection; numerical controls/output locations cannot authorize coarsening
or changing surface material/forcing. The actual pilot receipt must report
the matching method, TDB epoch, convergence metrics and complete native CSV.
No accepted criterion, model epoch, thermal depth grid or review is filled in
by the Python verifier.

## Native state and numerical campaign binding

The exact CSV header is
`logical_surface_element,depth_index,depth_m,kinetic_temperature_K`. Every
logical element has boundary index -1 at depth zero followed by all layer
centres indexed 0 through layer_count-1. Duplicate/missing cells or layers,
inconsistent depth grids and nonpositive/nonfinite temperature are errors.
This is a full model state, not a selected-site average or interpolation.

State bytes alone cannot authenticate physical inputs: a uniform grid can
remain identical after terrain or ephemeris changes. Every native initial-
state receipt must therefore include
`thermal_initialization_spice_kernel_count` and zero-based
`thermal_initialization_spice_kernel_<index>_file` from CSPICE's actual
`ktotal/kdata` furnished pool. The verifier retains the **ordered actual
absolute filenames**, including repeated loads, reads their hashes and checks
compatibility with the frozen pilot manifest. It does not substitute a
predicted installation directory or merely compare configured paths. Loaded
file order/bytes and configured LOLA IMG/LBL bytes must match the qualified
pilot in every case. Changed kernel/DEM content is ERROR even when the full
temperature CSV is unchanged. Missing actual-pool provenance or a file/hash
record is SKIPPED. Additional I20 observation-file entries are permitted;
their separate data/operator gates remain authoritative.

Every formal case must request `m2_verification_record_initial_thermal_state`
and retain the actual native artifact. Scoring verifies all 288 states and
canonical staged physics identities before numerical residuals/intervals are
evaluated. The formal import aggregates each state path/hash and I20 verifies
its actual cohort state belongs to that same I19 evidence list.

Numerical controls excluded from the physical projection still have their
own frozen recipe checks. Each actual input retains the original template's
cohort packet scale and source-integral tolerance, and its ordered surface
elements must equal the independently frozen native-site registry (including
shared elements). The verifier does not derive those choices from the
already-executed cases.

The separate source normalization must also be regenerated through an explicit
`production_input_tier=cluster` pilot. Its actual input, native manifest and
receipt are path/hash bound; both actual and template physics identity must
equal the reviewed matrix input. An internally consistent old development
pilot cannot supply a different mesh/source/thermal problem's normalization.

The measured cost forecast must explicitly include actual initialization
cost for the complete selected case count. The implementation introduces no
free cached-state loading and must not use an unexecuted periodic pilot or a
generic I17 PASS to waive that work. Long execution remains opt-in, forecast,
memory, disk and total-allocation budget gated.

Concretely, the forecast's `initialization_cost` object must contain
`included_in_allocated_core_hours=true`, `reexecuted_for_every_case=true`,
the canonical `qualification_sha256`, `native_state_sha256`,
`initialization_mode`, `case_count=288`, `measurement_mpi_ranks`, positive
`measured_per_case_wall_seconds`, `allocated_initialization_core_hours`, and
nonempty absolute-path `evidence_files` descriptors. The initialization
allocation must cover the measured per-case time times all cases/ranks, and
the existing total allocation must include it. A documented conservative
allowance is permitted; free or unmeasured cached initialization is not.

The formal manifest retains the original `native_evidence` family of exactly
288 real receipts and adds `native_source_evidence` for every validated case
file, including the complete initial state. `applied_initial_state` identifies
the qualified state, physical identity and all actual per-case hashes. I20
checks its finest-cohort state's membership in this source evidence rather
than accepting a different CSV or using a synthetic aggregate receipt.

Local tests use clearly labelled temporary **protocol/state-grid arithmetic
fixtures**. They do not qualify a real initial state, create observations or
approve any new physics assumption. Native actual-epoch qualification remains
missing until an investigator supplies the required real artifacts and review.

## Separate I20 thermal-history initialization

The annual thermal history is not implicitly assigned the cohort experiment's
initial epoch. Its reviewed polar contract must explicitly name an
`initial_state_qualification` path/SHA-256 descriptor, using the same reviewed
state schema with real native pilot artifacts. The selected pilot can be I20,
or an explicitly reused I19 proof when the actual I20 physics/mode/epoch/state
match. There is no automatic lookup/fallback to a convenient cohort proof.

The fresh I20 native exporter writes
`receipt.txt.I20_initial_thermal_state.csv` after its declared initialization
and before advancing to any query. Its independent receipt fields identify
the actual applied state, element/layer counts and TDB epoch. A legacy explicit
I20 spin-up may establish a later applied epoch than the original uniform
startup; that actual state must have its own qualification. The I20 verifier
uses that explicit state proof and unchanged native convergence metrics,
without declaring that an annual interval equals the cohort horizon or phase.
Missing own proof is SKIPPED; contradictory bytes/mode/epoch are ERROR.
