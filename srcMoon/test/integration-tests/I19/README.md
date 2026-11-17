# I19 — Surface inventory and global mass-conservation accounting

The ordinary `test.py` runner automatically executes the existing
`test_budget.py`, `../../test_m2_protocol.py` and the new
`../../test_m2_formal_campaign.py` provenance and
dispatch contracts before invoking the shared scientific runner. Therefore
the literal I19 entry in `run_all_test.py` covers these infrastructure checks
even when formal qualification is SKIPPED. Commands, source hashes and each
child's stdout/stderr remain under `I19/protocol_checks/`. A child failure is
reported as I19 infrastructure `ERROR`, with native execution not started.
A passing bookkeeping/metadata fixture does not qualify the production
conservation campaign, replace native receipts or open a predecessor gate.

## Applied initial-state qualification

[M2_INITIAL_STATE_CONTRACT.md](../../M2_INITIAL_STATE_CONTRACT.md) specifies the
required reviewed `actual_epoch_thermal_qualification` schema and all native
receipt/state/forecast fields. A generic READY proof about the requested epoch
cannot qualify a state the code did not apply. The diagnostic default remains
`thermal_initialization_mode=declared-uniform`; explicit `periodic-spinup`
applies the existing production spin-up before source/weight queries and
requires its cycle end to equal `simulation_start_utc`. No interval or initial
condition is selected by the test runner.

Formal staging preserves the original mesh, AMR, material, source, thermal,
forcing and clock settings. A canonical physical-input identity excludes only
the explicitly named numerical/output controls. It rejects silent coarsening
or other physical changes. Every formal case records the actual initial
boundary and **all** thermal layers in
`receipt.txt.I19_initial_thermal_state.csv` before cohort seeding, even when
restart/checkpoint output is disabled.

The independent scorer authenticates all 288 complete applied states before
computing numerical differences. Each must match the explicitly qualified
pilot state SHA-256, method, UTC/TDB epoch, full element/layer grid and physical
input/binary/source identities. A qualified uniform state is permitted only
when its review applies to that actual uniform state; periodic evidence is
never attached to an unexecuted initialization. These are model initial
conditions, **not independent observations or validation truth**.

The actual CSPICE furnished-pool count and ordered absolute file list are
native receipt fields, not a Python prediction from `spice_path`. Their bytes
and configured LOLA IMG/LBL hashes must agree with the reviewed pilot and
actual case manifests. A uniform state SHA can match while these physical
files differ; that is still ERROR. Missing native loaded-kernel provenance
keeps old artifacts SKIPPED. Unrelated extra observational file entries do
not change the initial-physics comparison or weaken their own data gates.

Missing exact qualification or an older native exporter leaves I19 SKIPPED.
Changed states, incomplete grids or contradictory method/epoch/source/physics
identity are ERROR. Actual numerical misses retain FAIL; no threshold is
relaxed. The measured total cost forecast must include re-executing the
qualified initialization for every case, including its full-state diagnostic
output. No free cached-state loading is implemented.

The aggregate retains exactly 288 real receipt descriptors in
`native_evidence`, plus `native_source_evidence` for all validated case files
and `applied_initial_state` with the actual per-case state hashes. I20 binds
its cohort artifacts and states to that same evidence. The existing
`test_m2_formal_campaign.py` suite now includes explicit uniform/periodic
protocol fixtures, complete-grid and changed-state refusals, all-288-before-
scoring control flow, and initialization-cost checks; the ordinary I19 runner
executes it automatically. No real state qualification is created by those
temporary infrastructure fixtures.

## Supporting development before M1 closure

The user-approved November 12 development-order exception permits a separate
supporting run while formal I19 remains predecessor-gated:

```bash
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

This builds the dedicated Na/Ar executable, enabling individual weights for
partial finite-reservoir packets and compiling chemistry OFF. It runs local
U08-U13, read-only package QA, supporting I14/I18 adapters and the supporting
I19 budget fixture, then restores canonical Moon/SPICE. The normal formal
I14-I20 runners are also invoked to preserve their genuine I13-blocked SKIPPED
results. Supporting PASSes have their own counts and never become formal
I19 or M2 PASS. An existing M2 build may be explicitly supplied with
`--m2-amps /absolute/path/to/amps-moon-m2-na-ar` and matching generated tree.

`MoonM2Budget` observes the same production source, surface and mover paths
as ordinary runs. It independently sums owned gas particles, replicated
surface stocks counted once, and local pending adsorption-minus-release
fluxes. Surface deposit is distinct from whole-particle gas adsorption;
uncredited adsorption is a physical surface sink. External injection equals
accepted emission minus surface release. Signed prescribed replenishment is
separate from gas injection. Gas, surface, combined atom and mass identities
are retained by species and independently recomputed with Decimal arithmetic.
The floating-point error bound includes the actual double surface-array
updates as well as long-double sums; the 0.5% source criterion remains intact.

The fixture synchronizes both species to the smaller initialized global step,
performs one initialization step, adds 64 controlled incoming particles and
one outward high-energy particle per species, and seeds 32 model-weight
equivalents into one logical surface element before the initial ledger row.
These are explicit verification initial conditions, never observations or
production equilibrium claims. Twenty ordinary steps exercise actual sources,
adsorption, release and escape. The predeclared I09 classifier uses
`v_inertial=v_LSO+Omega_LSO cross r_LSO` and specific two-body lunar energy:
positive-energy exits are lunar escape; bound finite-box exits and known
mover removals are numerical losses. This is not a claim about ultimate
three-body escape. Both classifier branches are checked explicitly; numerical
loss during the ordinary fixture prevents a supporting PASS.

The `supporting/I19/` directory contains `m2_budget.csv`, `budget.png/.eps`,
input, receipt, source/build/executable/SPICE manifest, stdout/stderr and
resource history. `analyze_budget.py` preserves all rows and requires actual
Na/Ar emission, adsorption, deposit and escape plus Ar release, preventing a
vacuous zero-event PASS. The signed source/sink equations are documented in
`test/M2_IMPLEMENTATION_PLAN.md`; complete formal I19 remains unresolved until
its predecessor and complete scope pass.

Purpose: Close synchronized Na and Ar particle and mass budgets across adsorption, desorption, inventory change, escape, physical sinks, and numerical loss.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I19 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: No external observation data. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The approved M2 budget case uses neutral Na and Ar with photolytic chemistry
compiled off. Consequently `converted_in` and `converted_out` are required to
be exactly zero, not estimated from a noisy ensemble. The dedicated Na+Ar
build helper verifies the generated `_PIC_PHOTOLYTIC_REACTIONS_MODE_` macro and
records its SHA-256 before it preserves the executable.

The budget table must separately report physical atoms and kilograms for
external injection, adsorption, desorption, initial/final exosphere,
initial/final surface inventory, physical escape, and numerical domain loss.
Adsorption and desorption are internal transfers and cancel from the global
closure equation; they remain explicit diagnostic columns. Deterministic
bookkeeping must close to machine precision, and the full weighted Monte Carlo
window must have absolute imbalance below 0.5% of its dominant source.

Formal I19 is currently `SKIPPED / NOT VALIDATED` because same-campaign I13
is unavailable. The supporting eight-rank Na/Ar fixture now closes its actual
gas/surface identities and exercises weighted finite release, incoming
adsorption and outward escape. This scoped result does not qualify every
mover, source, chemistry conversion, hybrid layout or time/weight resolution.
Those required paths and the complete formal window must be established
before formal I19 PASS; existing surface flux arrays alone remain insufficient.

An executed formal comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.

## Approved November 13 numerical extension

The original 20-step positive-coverage fixture remains intact. A separate
optional full `m2_verification_*` input group selects source-fed, source-free
cohort or source-normalization adapters. Their shared production transport,
exact signed accounting, source exclusions, independent dt/weight controls,
actual seeds, indexed physical clocks, passive tag semantics and source-only
quadrature are detailed in [M2_NUMERICAL_PROTOCOL.md](../../M2_NUMERICAL_PROTOCOL.md).
The initial stock is specified in physical atoms independently of the refined
packet weight. Cohort fields migrate in real particle records and take no
extra RNG draws. First departure and current occupancy are separate metrics;
returned atoms do not reset first residence.

The global development phase runs all 36 short cases automatically, retains
per-stage CSV/PNG/EPS and counts them as supporting checks. Two short seeds
cannot supply full formal confidence or long retention qualification. The
explicit cluster tier exports the complete 288-case plan, then requires actual
predecessor/data/site/normalization/statistical/cost proof before execution.
The original exact count/atom/kg and <0.5% source accounting requirements are
unchanged. Approved 1% equivalence and 0.1% long-run numerical-loss budgets
apply only to their separately declared numerical metrics. No observational
error or physical coefficient is inferred from these allocations. Source,
preflight and supporting PASS never become formal I19 PASS by themselves.

## Formal campaign assembly and reviewed execution contract

`m2_formal_campaign.py` connects the existing cluster executor and independent
scorers to the normal I19 result layout. It adds no physical law, source term,
observation or acceptance relaxation. `m2_formal_contract_file` in the `.in`
receipt selects an investigator-reviewed JSON dispatch contract. An absent,
unapproved or unavailable contract remains SKIPPED before any model-result
inspection. Ordinary/development execution never implicitly launches the long
design: a contract selecting execution additionally needs the explicit
`--allow-long-execution` CLI option (the global runner's corresponding option
must be supplied explicitly when using that entry).

There are two actions:

* `import` revalidates an already completed, authentic **same-campaign** long
  numerical design and reruns its independent science scorers. Import is not
  a copy of its old PASS label. The selected worker directory must contain
  `job_plan.json`, `preflight.json`, `cluster_summary.json` and all 288
  `cases/<key>/` directories. Supporting/short results cannot meet this
  contract. A missing case stays SKIPPED; an altered identity is ERROR.
* `execute` dispatches `run_m2_cluster.py --execute` under the existing
  all-process resource/deadline supervisor. The approved physical horizon,
  sixteen seeds, nine dt/weight cells and two physical cases remain unchanged.
  Missing data, thermal initialization, statistical justification, measured
  cost or predecessors prevent launch. This work package has not run that
  long qualification.

The JSON schema is `srcMoon-I19-formal-execution-v1`. Required top-level
fields are:

| Field | Meaning / required value |
|---|---|
| `schema`, `status`, `action` | Schema above, `APPROVED`, and `import` or `execute`; the runner never writes an approval. |
| `review` | `approved: true`, a nonempty `reviewer`, timezone-qualified `approved_utc`, and `frozen_before_execution: true`. The freeze asserts that the physical design/support evidence was fixed before its native execution; dispatch review does not retrospectively approve different physics or selections. |
| `qualification_root` | Actual current campaign containing I13–I18, not an unrelated directory of copied PASS labels. |
| `protocol_sha256` | Exact current approved `M2_I19_I20_PROTOCOL.json` hash. |
| `production_source_signature_sha256` | SHA-256 of the complete `production_source_identity()` dictionary serialized as UTF-8 JSON with sorted keys and separators `(',', ':')`. This covers dirty/untracked production code, not merely Git HEAD. |
| `input`, `executable` | Each exactly `{"path": "...", "sha256": "..."}`; binds the caller's actual `.in` and retained native binary. |
| `sites_file`, `forecast_file`, `normalizations_file`, `statistical_support_file` | The same exact path/hash records for the original cluster preflight. A path alone is insufficient. |
| `execution` | Explicit `mpi_ranks`, `workers`, `allocated_cores`, `core_hours`, `memory_limit_bytes`, `artifact_limit_bytes`, `host_free_floor_bytes`, `disk_free_floor_bytes`. No resource estimate or statistical claim is invented. |
| `completed_summary` | Import only: exact path/hash record of the completed **worker's** `cluster_summary.json`; its sibling `cases/` must contain the full native design. |

Paths inside the JSON resolve relative to its own directory. The `.in` path
resolves relative to the `.in` file. JSON duplicate keys/non-finite numbers,
changed hashes, unsupported actions, contradictory MPI/resource settings and
stale native-source signatures are ERROR. The builder verifies the real
generated chemistry/RNG layout and native binary before qualification. The
existing cluster preflight still controls all required scientific/data/site
and cost prerequisites; this adapter does not bypass them.

Formal input staging preserves the reviewed template's AMR, source switches
and strengths, thermal material/time/initialization controls, and epoch. The
I18 bounded-volume/Ar-enabling transformer is used only by the development
tier. Formal input must already declare Ar enabled, an evolving environment
and zero initial surface inventory; otherwise the runner reports SKIPPED.
Canonical `physics_identity_sha256` binds the actual staged physical settings
to the reviewed template while allowing only named seed/numerical controls
and output suppression/ownership to differ. Unknown settings are retained in
identity rather than silently ignored by a prefix rule. Relative input paths
retain their assigning-file meaning, including `#include` fragments.

Each cluster case also requests the actual native initial full-column CSV and
records a prelaunch UTC timestamp. Qualification must authenticate the state
the native process installed before cohort/source injection. A separate
periodic spin-up report that was never applied does not qualify that run.
Omitted initialization mode retains the production `declared-uniform` default
and does not itself establish a physically qualified lunar initial state.
The detailed staging contract and parser-only verification are in
[FORMAL_STAGING_CONTRACT.md](../../FORMAL_STAGING_CONTRACT.md) and
`srcMoon/test/test_m2_protocol.py`; no observation is used by those fixtures.

Prepare a formal deterministic model-input normalization with
`prepare_m2_normalization.py --tier cluster` (the complete command is in
[M2_NUMERICAL_PROTOCOL.md](../../M2_NUMERICAL_PROTOCOL.md)). The default
development tier intentionally retains the bounded fixture and cannot be
promoted to formal use through an original-template hash. READY references
retain the actual staged physics hash, actual input and native manifest/
receipt descriptors, plus the initial per-site stocks and cohort-packet count
from which the ordinary-initialization reference weights are derived. Thus
the reference cohort stock uses the same AMR-dependent native weights as the
reviewed transport problem. Missing new provenance remains SKIPPED until an
appropriate reference is generated; no long run or reference approval is
performed automatically.

The normalizer's `--record-initial-state` option exercises the actual I19
initial-state writer before its source-only early return. It verifies the
full native grid and actual method/epoch/SPICE/DEM provenance, then retains
an explicitly unqualified candidate. This supplies code-path evidence and
review material, not initial-state qualification or observational evidence.
The global development runner requests it during source normalization;
requested missing or corrupt output is ERROR rather than an omitted artifact.

For every imported case the adapter checks actual native seed/layout,
executable, full physical duration, independent dt and weight controls,
configuration snapshot and external-data hashes. It independently recomputes
signed atom/kg conservation and cohort hazard/tag identities before applying
the existing sixteen-seed numerical-equivalence scorer. Matching resolution
means alone would not establish correct conservation, so both layers are
mandatory. Retained FAIL cases are never silently promoted by a new summary.

Once a real reviewed contract and all predecessors exist, direct dispatch is:

```bash
python3 -B srcMoon/test/m2_formal_campaign.py \
  --amps /absolute/path/to/verified/amps-moon-m2-na-ar \
  --input /absolute/path/to/reviewed/qualification.in \
  --qualification-root /absolute/path/to/current-campaign --mpi-ranks 8
```

Add `--allow-long-execution` only for an explicitly reviewed `action=execute`
contract and provisioned allocation. No example READY contract is shipped:
its physical inputs, source/build identities, approved evidence and measured
costs must describe actual retained files rather than invented placeholders.

Canonical evidence is written under `<qualification-root>/I19/`: input and
execution-contract snapshots; stdout/stderr; `formal_preflight.json`;
`run_manifest.json`; `result.json`; `native_evidence_manifest.json`;
`aggregate_native_receipt.json`; independent numerical score and evaluated
`comparison.csv/.png/.eps`. The numerical figure compares normalized
differences with zero, the analytical equivalence null; error bars are the
approved simultaneous **numerical confidence** intervals, not experimental
uncertainties. A skipped comparison produces no figure. No artificial
single native `receipt.txt` is created: the manifest's `native_evidence`
list binds all 288 real per-case native receipts by path and SHA-256.
Existing meaningful result directories cannot be overwritten.

### Control-plane verification

```bash
python3 -B srcMoon/test/test_m2_formal_campaign.py
```

The tests cover missing contract/predecessors, stale input/source signatures,
explicit long-job opt-in, same-campaign identity, layouts/resources, duplicate
or non-finite JSON, changed seeds, incomplete case families, false cached PASS
summaries, immutable prior results and exact arithmetic-row serialization.
All fixtures are temporary metadata/arithmetic controls, never observations,
native model output or evidence of physical qualification. They do not launch
MPI, inspect December residuals or alter canonical data. Their PASS is
infrastructure verification only; formal I19 remains SKIPPED until the real
long campaign and original scientific gates pass.
