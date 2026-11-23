# Stage 11 common material-background contract and test plan

The shared implementation is
[`sep_material_background.h`](../../src/models/sep_common/sep_material_background.h).
It is header-only and preserves the seven-member canonical archive and existing
BackgroundSnapshot layout/symbols. The
[stage-10 audit](STAGE10_BACKEND_CAPABILITY_AUDIT.md) freezes the current profiles
and missing physics; none is promoted by this interface work.

## Delivery status (2026-11-23)

**Stage 11 is complete for the approved scope.** Both maintained Make
qualification targets pass. Stages
10–11 retain the existing reduced SSE/Parker exterior as the only qualified
material profile. The exact twelve-profile inventory, immutable owner/interval,
canonical SI exports, pre-publication binding and fixed-step admission are
implemented; no missing model physics is supplied by a capability declaration.

Evidence root: `/data/vtenishe/TEST/sep-stages10-11-20261123T030052Z`.
The pre-integration matching builds and campaigns use the `region` names below.
The final integrated builds and fresh validation use the `integrated` names in
the follow-up record. Earlier
`*-gen0` campaigns completed 48 native, 58 material, 44 checkpoint and three
400-step applications plus eleven audit controls successfully. Those results
predate the final trajectory-region correction and remain historical evidence.
They are not substituted for current-source qualification.

Three real defects were corrected without changing the existing physical law,
RK controls, fixed AMPS timestep, independent reference or acceptance bounds:

- Unpopulated mandatory vectors/tensors could appear as plausible zeros. They
  now remain explicitly invalid until populated; MBC08 rejects them. The
  `uninitialized-valid-before/after` witness retains the actual failed/passing
  executions. These construction guards never replace missing physical data.
- Initial registered PIC geometry legitimately has generation zero. Explicit
  generation presence distinguishes it from missing metadata, while native
  installation independently checks actual registered ownership. MBC03/MBC05
  and the native pre-registration rejection test cover both requirements.
- The new trajectory wrapper accepted PFSS coronal points although its state
  API correctly rejected them as outside ParkerExterior. The wrapper now
  validates owner identity and actual region at the initial point, every
  existing RK evaluation and the final accepted point. This also covers zero
  duration and reverse-time queries without assuming monotonic trajectories.
  Legacy MoveVertex is unchanged. MBC18 reproduces 19 PASS/1 FAIL before the
  correction and 20 PASS afterward, retains bitwise-positive trajectory
  equality, and independently achieves relative error 2.72e-16 against the
  unchanged 2e-12 reference bound.

Executed results before the fixture integration:

| Matching evidence | Actual result |
| --- | --- |
| `contract-region-strict`, `contract-region-sanitizers` | 20 PASS each; strict warnings and ASan/UBSan, actual exits 0 |
| `shared-material-region` | All 20 cases PASS through the shared regression registry |
| `pcli-region` | PCLI01–12 PASS on two MPI ranks with two repeats; actual exit 0 |
| `native-region-tests` | All 48 jobs exit 0; all 15 IDs PASS at ranks 1/2/4/8; collector exit 0 |
| `backend-native-region` | SBACK06–08 PASS on two ranks, actual expected exit 1 before allocation |
| `controlled-region/results` | CV01–12 and IV01–06: all 18 PASS, actual launcher exit 0 |
| `material-region` | All 58 restart/repartition roles PASS, collector exit 0 |
| `checkpoint-region` | All 44 roles PASS; byte-exact continuation up to 21,425,968 bytes; collector exit 0 |
| `application400-region` | All three full applications and original independent audits PASS; eleven corruption controls PASS |
| `bounded-region` | Full Make exit 2: bare native --all-tests cannot supply heterogeneous native prerequisites; complete 222-ERROR catalog retained |
| `region-harness-compatibility.json` | Orchestrator 51 PASS, Python runner 27 PASS, catalog/traceability PASS, retained public ABI/header and seven-member archive checks |

The final `material-region` campaign passes all 58 roles, and
`checkpoint-region` passes all 44 roles with a maximum budget residual of
1.466477706182233e-15 against 1e-10. All three `application400-region`
applications exit 0 and pass the independent verifier: each completes 400
unchanged 300-s steps, all 357 surviving segment exposures reconstruct exactly
(5–400 samples), and the maximum relative budget residual is
2.2022684795608913e-15 against 1e-10. All eleven original corruption controls
also pass. Original 1800-s application and 300-s checkpoint deadlines remain
unchanged. No preparation or partial run is used for these verdicts.

The full `bounded-region` Make invocation nevertheless exits 2. Its source
prerequisites and controlled comparisons pass, but the bare native --all-tests
recipe chooses the maximum initialization requirement for a heterogeneous
catalog. The default legacy lines have no registered material owner, so the
unchanged production-step prerequisite rejects before callbacks and records
222 ERRORs. This is an actual remaining gate failure, not a passing routine
suite or 222 independently evaluated physics failures.

The executable's exact routine inventory is 193 initialization-free cases and
29 native field-line cases. A diagnostic invocation selecting exactly the 193
initialization-free IDs exits 0: 175 PASS and 18 explicit CV/IV SKIP without
native manifests; those 18 separately pass the approved full-reference workflow.
This diagnostic does not qualify or remove the other 29 cases.

A descriptor-aware isolated routine driver and explicit native fixtures were
developed in an evidence-owned prototype. The 16 native storage/accessor
cases need allocation rather than a production particle-step preflight; five
CP05C cases can reuse the qualified reduced-material protocol-2 setup. FTE01
and PARKER01 declare zero centerline flow and need actual scoped identity-map
ownership after installing their manufactured fields/tensors. BIRTHN01–06 need
fixture/lifecycle migration, retaining their original duration and byte-image
oracles and an actual compatible checkpoint for the negative save test.

BIRTHN01–05 are explicitly operator-stage software tests. The maintained method
freezes departure geometry/coefficients during particle transport, then advances
the material geometry after readers close. Their repair preserves that existing
scheme, the declared nonzero radial wind, exact 0.0025-s remaining-duration
comparison, 0.01-s parent/full-host negative control, and the large-epoch test's
explicit relative durations. It must not add a private dense-time physics branch
to the movers, silently restrict curvature, mutate a live AMPS timestep, publish
a false material owner, or classify missing fixtures as SKIP.

Direct `mpiexec -n N --bind-to none` launches work. The earlier general claim
that MPI could not run was incorrect: it arose from wrapped startup failures.
The retained stage-8/9 SWMF receipts already show nine successful MPI phases and
nine original-reference numerical comparisons. Rank logs now use MPI's own
`--output-filename`. A sequential validation driver needs independent child MPI
sessions: simply reusing one parent rank fails after the first child finalizes.
The evidence-only `independent_runner.py` clears inherited MPI session metadata
before the unchanged driver, and preserves site tools on PATH while selecting
the installed scientific Python. This local execution fixture changes no model.

The four `native*-region` builds are separate copies preserving the older build
evidence. Their physical PIC/thread macros remain byte-exact. Seventy actual
PIC-consuming objects in each copy were recompiled to regenerate dependencies
against that copy's local PIC headers; dependency files were never rewritten
to manufacture provenance. The unchanged maintained audits accepted the actual
compiler/object/header/archive closure before native execution.

The full-source run also exposed stale validation plumbing. CV01 now hashes the
canonical relocated transport implementation. The relocation gate requires all
six canonical include roots exactly in both supported layouts. Two mode-0644
POSIX scripts use Make's shell; runner tests retain the exact script/target
checks. SCOPE02 admits solely the Cartesian cache derived from the committed
field-line coordinate, with fifteen rejecting mutations and native SINT08 still
required. The turbulence audit requires exactly one live call in the library's
actual amps_time_step and none in the driver, with four rejecting mutations.
Pre-existing ignored bytecode was preserved outside the source tree; DOC03's
clean-source requirement remains intact. These are faithful contract repairs,
not removed numerical assertions or relaxed ownership checks.

The broader shared SEP-common baseline records 91 PASS/52 explicit SKIP.
ASan/UBSan runs use detect_leaks=0 where LeakSanitizer cannot operate under the
traced session; no leak-coverage claim is made. Optional STATIC01 reports SKIP
when this compiler lacks -fanalyzer. No observational validation is claimed.
The stage-10 full-regional original 180-s MASS-TIME timeout remains ERROR and
that profile stays unsupported. A separate pre-existing packaging CLI mismatch
also remains ERROR: the Make target requests --root --self-test, while the
unchanged overlay helper accepts a positional root/--allow-generated. It is not
included in a passing packaging claim or silently bypassed.

## Maintained regression integration and matching builds

The final prototype passed all 222 Routine cases and all 24 native Birth roles
over two seeded repetitions. The reviewed C++ changes are integrated. The four
new configured builds recompile their PIC consumers and match the applied source
hashes. The maintained Make targets now use separate startup fixtures for the
different native prerequisites. The exhaustive runner applies the same routine
fixtures and keeps every Extended descriptor and its result. The independent
controlled comparisons retain their original input decks and reference checks.

The focused and application gates on the new builds passed as follows.

| Evidence under the same root | Actual result |
| --- | --- |
| `contract-integrated-strict`, `contract-integrated-sanitizers`, `shared-material-integrated` | All 20 material-contract cases pass in each campaign |
| `pcli-integrated`, `backend-native-integrated` | All 12 production CLI cases and three intended two-rank unsupported-backend rejections pass |
| `native-integrated-tests` | 48 successful MPI processes cover all 15 required IDs on ranks 1, 2, 4 and 8 |
| `material-integrated`, `checkpoint-integrated` | All 58 material and 44 checkpoint/threaded roles pass, including byte-exact continuation up to 21,425,968 bytes |
| `birth-repeat-integrated`, `orchestrator-integrated`, `selected-cp05-integrated` | All 24 repeated Birth role jobs, the ten real orchestrator cases, and the selected CP05C01 Make workflow pass |
| `application400-integrated` | All three movers complete 400 steps of 300 s and pass their independent verifiers and eleven corruption controls |
| `controlled-integrated` | `make test-controlled-analytical` exits 0 with CV01–12 and IV01–06 all PASS |
| `bounded-integrated` | `make test` exits 0 with all 222 SEP Routine cases and four embedded SWCME Routine cases PASS, plus every source prerequisite |

Each full application reaches 120,000 s with 357 sampled segments and exposures
from 5 to 400. The largest independently verified relative budget residual is
2.203e-15 against the original 1e-10 bound. Two checkpoint launch receipts could
not be saved on their first attempt because the evidence capture command
exceeded the system argument limit. The raw logs were preserved and the same
two roles reran successfully under the original deadlines. The final collector
uses those actual rerun receipts.

The final audit `stage10-11-qualification.json` checks exact Routine membership,
all 222 actual case exits, matching executable and fixture hashes, and all
current build and campaign summaries. The main executable SHA256 is
`11c17e4be76b894f843b7a1152e74b26b97dd3bbb36a678d701d0daf0e7df138`.

## Identity, ownership and availability

The descriptor declares exact profile/revision; physical event/configuration;
optional native run binding; separate background, geometry and native snapshot
generations; material lineage; composition/EOS; HCI/SI; density basis and closed
physical coverage. A cache generation is not a material label or native
snapshot generation. Proper density, frame density and core nonrelativistic
fluid density are distinct declarations; the interface performs no guessed
conversion.

Geometry-generation presence is explicit: zero is the valid initial generation
of a registered PIC line. Absence rejects a material trajectory. Registration
identity is checked independently in the native adapter, so admitting zero
does not admit an unowned geometry or alter PIC's sampling-generation counter.

The owner's maximum possible capability mask is bounded by its audited profile.
An actual owner declares a subset, and each query states required quantities.
Primitive, species, thermal, full grad(U)/grad(B), continuous-time trajectory,
inverse/deformation, shock/contact and downstream-volume capabilities remain
separate. Absent optional fields are explicit. A present NaN or invalid value
always rejects; zero alpha abundance is a valid declared population, not missing
composition. EOS semantics remain owned by the provider; canonical export
checks independently verify the current ideal-gas composition.

Mandatory vectors/tensors are invalid until every component is populated.
Their NaN construction guards are rejected, never interpreted as physical
absence or replaced with zeros. Optional physics still uses explicit optional
records. MBC08 covers forgotten U/B/tensor population; the isolated
`uninitialized-valid-before` witness fails before this correction and
`uninitialized-valid-after` passes with the same omitted-data input.

An immutable Interval retains its actual const owner and captured descriptor.
Queries reject changed ownership metadata before/after callbacks, wrong epoch,
position, side/region or generation, unavailable data, unknown profile/revision
and out-of-coverage times. Full physical-step admission requires exactly
end = begin + configured fixed dt and complete coverage; it never clips a
step to the available interval. Branch-preserving derivatives are required;
a stored divergence or density history cannot replace all nine tensor entries.

The reduced adapter exports actual canonical SI samples and the existing
endpoint-preserving trajectory method. A captured owner survives later
reconfiguration without adopting new physics. Surface-local RH limits are
explicit one-sided primitives; missing downstream temperatures/gradients and
a downstream volume stay unavailable. No new map, heat law or inverse is
invented for another backend.

Native binding validates the existing SnapshotStore provider/ownership,
native background configuration identity, coverage and explicit native
generation before its original
PreparePublication/CommitPublication transaction. The store remains the sole
native writer/read-phase authority. The reduced MaterialView owns a derived
const interval; checkpoint restore reconstructs it from the real saved owner.
It is not an added opaque pointer or physical payload in schema-5 images.
The native identity is the existing Background::CurrentConfigurationFingerprint,
not a replacement for the full runtime/checkpoint configuration hash. Metadata
reuse validates the actual retained store object, preventing a hypothetical
replacement from concealing a foreign live configuration. MBC20 retains that
negative control and verifies that rejection does not publish anything.

## Serialization and compatibility

The new SEP_MATERIAL_DESCRIPTOR version-1 envelope serializes identity and
capabilities only, with classic locale, binary64 round-trip precision and the
existing deterministic identity checksum. Unknown schema/revision/capability,
truncation, changed digest and trailing data reject transactionally. This is
not a physical-state checkpoint and does not authenticate external data.
Decoding it cannot create a missing physical owner or confer admission.

The existing two-argument bundle LoadAndPublish ABI, schema-1 neutral codec,
SWCME public records, legacy backend-inventory overload and BackgroundSnapshot
ABI remain intact. Production SWCME checks use the exact region mode selected
by the canonical configuration; no family-wide success is substituted.

## Focused gates

MBC01–20 are new software/adapter test identities; they do not replace
SBACK01–08, SINT, checkpoint or full-application gates.

| ID | Independent assertion |
| --- | --- |
| MBC01 | Exact twelve-profile catalog and stable names/IDs; one admitted material profile; unknown profile rejects |
| MBC02 | Actual canonical SHOCK_ONLY/FULL_ICME selection; legacy/planned, unsupported region/use/revision rejection |
| MBC03 | Exact descriptor round-trip including valid initial geometry generation zero and explicit presence; decoded metadata cannot construct a missing physical owner |
| MBC04 | Corrupt, truncated, foreign-version, trailing and forged-capability descriptors reject |
| MBC05 | Missing frame/units/basis/composition/EOS/lineage/generations reject before capture |
| MBC06 | Closed interval endpoints work; next representable outside time rejects before callback |
| MBC07 | Missing full gradient, species and populated temperature stay unavailable |
| MBC08 | Invalid present optional primitives/species/tensors reject even when not requested |
| MBC09 | Foreign event/configuration, time, position, region/side and generations reject |
| MBC10 | All nine independently specified tensor components survive; one-sided queries require identity; unknown side/region rejects |
| MBC11 | Complete fixed step accepted; partial, oversized, uncovered and absent-dt requests reject |
| MBC12 | Distinct generation roles bind correctly; foreign provider/read-phase publication leaves the actual store unchanged |
| MBC13 | Mutated owner metadata cannot rebind an already captured interval |
| MBC14 | Missing inverse/contact/front/downstream capability rejects before physical callback |
| MBC15 | Actual neutral export equals independently constructed canonical primitives and full derivative tensors |
| MBC16 | Original physical owner/data remain intact after a different-composition reconfiguration |
| MBC17 | Independent charge neutrality, ion mass sum and species ideal-gas pressure detect electron/proton-density confusion |
| MBC18 | Exact positive legacy trajectory and independent accelerating-wind flight-time oracle; actual PFSS, zero-time, reverse-crossing and foreign-owner requests reject |
| MBC19 | Actual front speed/RH limits retained; local jump never admits a downstream volume or invents thermal data |
| MBC20 | Actual owner/native binding and descriptor restore; late foreign realization leaves the store and valid owner intact |

Manufactured contract fixtures are explicitly labelled software tests.
Composition variations are declared verification inputs copied to isolated
directories; no production event, observation or reference is rewritten.
Source builds use -Werror and separate ASan/UBSan checks, with LeakSanitizer
disabled where unavailable. They do not replace native qualification.

## Native and routine completion gates

Rebuild the ordinary reduced field-line application with complete dependency/
object/archive identity. Run PCLI01–12 before/after, actual two-rank SBACK
preflight, native geometry/source/three-mover checks, material publication and
fresh restart/repartition/thread cases relevant to the new derived binding.
The final three-mover application must complete the unchanged 400 x 300-s
window with original independent verifier and eleven corrupted-output controls.
No source-only interface PASS substitutes for those executions.

Register this suite in Make/source selection, the complete SEP orchestrator
with a strict full-catalog reader, and the shared extended core registry.
Synchronize the generated catalog, traceability and pinned harness membership.
Missing compile/output/case evidence is ERROR with the full catalog and a
nonzero process result; prior evidence directories must remain untouched.


## Executed commands and retained results

Use fresh output directories. Every prepared plan retains its exact canonical
argv/cwd and source/build hashes; actual process receipts additionally record
the launch argv, true exit code and rank logs. The current evidence root is
stated above.

```sh
EVIDENCE_ROOT=/data/vtenishe/TEST/sep-stages10-11-20261123T030052Z
python3 -B srcSEP/test/run_material_background_tests.py --output-dir "$EVIDENCE_ROOT/contract-region-strict"
python3 -B srcSEP/test/run_material_background_tests.py --sanitizers --output-dir "$EVIDENCE_ROOT/contract-region-sanitizers"
python3 -B srcSEP/test/run_native_standalone_integration.py --build-root "$EVIDENCE_ROOT/native-region" --output-dir "$EVIDENCE_ROOT/native-region-tests" --ranks 1 2 4 8 --collect
python3 -B srcSEP/test/run_material_checkpoint_tests.py --output "$EVIDENCE_ROOT/material-region" --collect
python3 -B srcSEP/test/run_stage7_checkpoint_tests.py --output "$EVIDENCE_ROOT/checkpoint-region" --collect
python3 -B srcSEP/test/run_standalone_application.py --output-dir "$EVIDENCE_ROOT/application400-region" --collect
python3 -B srcSEP/test/test_standalone_application_verifier.py --campaign-dir "$EVIDENCE_ROOT/application400-region"
python3 -B src/test/generate_test_catalog.py --check
python3 -B src/test/validate_traceability.py
```

The commands containing --collect are completion gates, not assertions that an
in-progress campaign has finished. Exact current results are recorded in the
table above and each campaign's completed summary. Final Make/controlled launch
plans and receipts belong to bounded-region and controlled-region respectively.
All successful MBC records include actual assertion counts, wall duration,
seed, configuration/source provenance, exact metadata tolerance and the
unchanged SI-reference bound/residual where a numerical comparison occurs.
Cases without a numerical comparison do not invent a zero residual.

Older final-evidence-audit.json and the failed wrapped-launch/source-gate
attempts remain retained history. They predate the corrected direct-MPI method,
geometry-zero handling and final trajectory-region guard; none substitutes for
current matching-source execution.
