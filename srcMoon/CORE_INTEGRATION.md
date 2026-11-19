# Integration with the current AMPS core

## November 19 integration and corrective verification

Master integrates `60a156ddcc3deeaf48e2981b0beb14738245457e` with the preserved local source (eleven incoming commits). The completed corrective CORE registry has **844 PASS, 240 SKIP, 1 FAIL and 1 ERROR**. Exact historical/current-upstream PASS preservation audits pass; the known nonpasses remain visible. Moon M0 is **10/10 PASS**, with a zero-byte reviewed nightly diff, and **1533 targeted native checks PASS**. The U-series has 15 PASS/11 SKIPPED and zero FAIL/ERROR.

The original concurrent CORE attempt is **ERROR**: two CampaignBudget calls in
threads shared the same parent PID and adopted each other's children. Its
conservative outer reserve interrupted the merged lane before aggregate JSON
publication. Completed owner logs do not supply invented individual PASS rows.
The original failed deadline, logs and supporting short-tail retry are retained.
The completed pristine baseline is authentic: 843 PASS/240 SKIP/1 FAIL/1 ERROR.

A separately recorded corrective campaign uses one OS supervisor and one full
CORE candidate, with 24 build jobs on 36 logical cores. It retains its own ORIGINAL
3600-second deadline including setup/build/run/report/delivery, without resetting
that clock after failures. No owner timeout, fixture, reference or acceptance
threshold changes. Exact case-identity audits protect both historical 844 PASS
and current-upstream 843 PASS, and their union; incomplete/preflight evidence never
substitutes for actual runtime rows. Known incomplete-coverage FAIL and full-CME
180-second ERROR remain separate from preservation PASS.

The pristine CORE baseline and complete Moon runs retain actual 822dec07 source
revisions. Late upstream 60a156dd changes only an Earth source-contract Python
test; five focused checks pass. All native/CORE/Moon inputs and registry bytes
are unchanged by that addition. The corrective candidate executes on 60a156dd.
The 742 native-input hashes and the late-delta proof preserve these distinct bases.

Three-way merges retain upstream APIs and Moon's existing opt-in force-tracer
declarations in `src/pic/pic.h`, plus the local registration in
`MakefileTest/Table`. Its build uses `$(TESTCOMPILEJOBS)` as required by upstream.
The nightly planner now honors an explicit Table over an implicit cached
Makefile manifest, retaining explicit configured-tree priority and host checks.
Four new precedence controls, all 83 nightly unit tests, the official registry
check and six field-line helper cases pass. Three existing maintained chemistry
build helpers are explicitly unignored so future commits include required source.

Builds use fresh generated trees and normal Config.pl/GCC 9.4/MPI/SPICE. Native
manifests retain exact commands, input/configuration/header/executable hashes,
consumed data hashes, seeds, rank layout and stdout/stderr. Targeted Moon uses
one/eight MPI ranks; the reviewed nightly uses four MPI ranks/OMP 1 and produces the exact
reviewed data hash. I04's original generic scorer manifest is retained alongside
a SHA-bound actual-runtime supplement, rather than relabelling an unexecuted
command as the actual invocation. No observed random seed is invented.

Fast-forward delivery preserves unrelated local source, new files and all eight
site configuration files. No staging, user commit, push, reference replacement
or raw/installed-data edit occurs. Active generated build/executable and site
Makefile remain unchanged; remove build and configure/build normally before
using an active application executable after source delivery. Original source,
index and overlap backups remain available with hashes and recovery notes.

Formal M1/M2/M3 qualification remains open. Passing software fixtures and ideal
collectors supply no missing physical inputs or observational agreement. The
source/data roadmaps retain their different recorded hashes; approved source
resource/development amendments govern this work while original physical/data,
calibration/hold-out and predecessor gates remain binding.

Final evidence:
`test_output/srcMoon/upstream-core-merge-corrective-20261119T025503Z/completion_report.json` and `.md`.
The original error/native package is
`test_output/srcMoon/upstream-core-merge-20261119T015550Z/`;
its `latest-target/` holds late-delta and actual-application evidence. Actual
commands and hashes are in the individual native run manifests and package
contracts; preserve their source bases when reusing evidence.

## November 17 integration

The user authorized integration of upstream
`4bff09213b1bfa3f59e113ba8cbb2d56f22d30ef` with local changes based on
`17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`, preferring upstream code and
requiring preservation of passing core tests. The final target advanced to
`5ce92cb2feecadd11a1d60ace69153a7c7da4122` during testing. Its two extra
commits affect unrelated comet Python/tooling and the Earth nightly command;
all 656 native/config source identities match the tested 4bff092 basis. The
latest nightly unit suite and official check separately PASS. Raw physics
receipts retain their actual 4bff092 revision. No new commit, push, reference
replacement, acceptance relaxation or new physics choice is part of this work.
The durable evidence directory is
`test_output/srcMoon/upstream-core-integration-20261117T185033Z/`.
Its `completion_report.json`/`.md`, delivery journal and source snapshots
record the exact delivered state; [STATUS.md](STATUS.md) records completion.

### Upstream code retained

All upstream changes to cut-cell geometry/measures, particle movers and lists,
runtime generic surface APIs, sampling and macro insertion are retained.
The existing five dirty shared C++ files add optional PSD geometry and passive
finite-release/cohort observers; no further shared C++ fix was needed. Earlier
committed Moon finite-inventory, moving-terrain and event hooks are already
present in upstream. The new force-free surface mover cannot replace Moon's
acceleration-aware moving-terrain route without losing required physics.
The detailed source/caller/switch/test audit is retained with the evidence.

Two integration issues were corrected in application/build adapters:

* The M2 derived input switched `PhotolyticReactions` off while retaining an
  explicit application macro that switched it on. Upstream now correctly
  inserts application settings inside the guard, after the earlier setting.
  The old gate found an early `OFF` string and falsely accepted that build;
  the linked I19 guard refused it. The Moon builder now changes both controls
  for the declared chemistry-free M2 fixture. Canonical Na chemistry stays on.
  It preprocesses the complete `pic.h`/`Moon.h` include chain with the actual
  make-expanded compiler/flags before and after compilation, requires all
  relevant macros, and binds effective modes, four immutable header snapshots
  and recipe identity to the executable. Old incomplete/regex-only receipts
  are rejected. The original ERROR remains retained with the corrective run.
* Enceladus's exporter included absolute checkout paths in three provenance
  comments. Both pristine upstream and the first candidate reproduced the
  exact-byte U00 failure in an isolated checkout. The exporter now writes
  repository-relative paths for repository-owned inputs, retaining explicit
  absolute external paths. Both examples were regenerated by their owner;
  all science/hash/other bytes are unchanged. The exact comparison remains,
  and relocation/external-path controls are included in the global runner.

No observation, frozen physical coefficient, numerical tolerance or I04
reference was changed. Future retained M2/M3 reuse intentionally requires the
tightened build identity; historical receipts are never backfilled.

### Regression scope and configuration

Pristine upstream and the integrated candidate each completed the same 1085
CORE cases: **843 PASS, 240 SKIP, one existing strict-coverage FAIL and one
existing CME 180-second ERROR**. Registry-aware comparisons against both the
historical and current upstream reports preserve every PASS and selected ID,
with no new failure or status change. This is preservation, not an all-green
registry. The older 919-PASS aggregate additionally included 62 Enceladus
application cases and fourteen Moon U aliases. All **62/62** Enceladus passes
and **14/14** Moon aliases are preserved, completing **919/919** historical
PASS coverage. E3 completed 10,000 steps with zero H2O/CO2 ledger residuals.
Their separate evidence is reconciled in the completion report rather than
silently dropped. This does not assert that every historical nonpassing
application/convergence case was rerun or qualified.

Five additional upstream native mesh/surface targets PASS (34 reports with
643 PASS report lines, including summary lines); 64 nightly-infrastructure
unit tests PASS. The I01-I10 baseline PASSes with unchanged reviewed reference
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`
and zero-byte `test_Moon.diff`. One-/eight-rank actual Na event budgets and
nine-window storage checks PASS. Supporting M2/M3 results retain their limited
scope: formal M1/I13 and M2/M3 qualification remain open.

Builds and runtime work belong to isolated `/tmp/amps-upstream-core-integration-`
checkouts. The active application/configuration and old generated executable
are not rebuilt as a side effect of source delivery. The fast-forward temporarily restores only the backed-up
dirty tracked paths, then reapplies their reviewed merged bytes. Untracked
source is never cleaned; original source and index backups are retained.
Start the next active
application build by removing `build/`, then use normal Config.pl and make -j.
The tracked Makefile receives upstream's input prerequisite and retains the
local SHARE linker option; it remains uncommitted. Generated `Makefile.test`
receives upstream registrations/recipes and retains eleven local fast-wave
include overrides. Other site configuration bytes are preserved.

Upstream's small no-receipt nightly mesh output and final particle-only
`ParticleData.restart` are retained for legacy regressions. Configured science
mesh output still follows its input, and science/I11 autosaves remain off.
Neither diagnostic output nor the legacy particle restart qualifies complete
surface/thermal/RNG/sampler restart or any observational comparison.

## Historical November 12 integration

The user authorized this work on 2026-11-12 without Git commits or pushes.
The pinned upstream revision is
`09fabea55e367f07a530460b488fa7a981ff6176`; the previous Moon development
checkout was `22ad748432ce032775189ba234be7a3d339a4d8d` plus preserved local
source. The tested source has been delivered into
`/data/vtenishe/Moon-test/AMPS`: master fast-forwarded to the existing pinned
upstream commit and all Moon/integration changes remain unstaged. No commit,
stash, reset, blanket clean or push was used. Original source/configuration
backups and intermediate evidence remain outside the code tree at
`/tmp/amps-moon-core-integration-20261112T224309Z-ehx53mh8/`.

## Integration contract and production choices

The upstream Poisson waiting-time, cumulative source selector and surface-angle
samplers remain the production authority. Moon's finite-reservoir release is
kept separate from external stationary injection, so no released inventory is
injected twice. No legacy rejection sampler is restored. Upstream generic
container, mesh, RNG, field, transport and performance work is preserved.

Necessary shared changes retain the production sampling-window reset, optional
actual-particle event observers, optional runtime source/reservoir parameters,
moving departure epoch/force-mover route and finite-release injector. Their
default/null behavior remains the normal upstream mode. A NaN departure scope
does no shared writes/reads, avoiding a new race in non-moving parallel sources;
finite moving injection keeps its explicit threading restrictions. The finite
release physical law, weighted packet representation and existing Moon
scientific assumptions remain those documented in the M2 contracts.

| Shared production location | Retained change and reason |
|---|---|
| `src/models/exosphere/Exosphere_SurfaceSource.h` | Unchanged upstream samplers: endpoint-safe Poisson waiting times, cumulative source selection, corrected angular distribution. |
| `src/models/exosphere/Exosphere.cpp` / `.h` | Optional finite-inventory dispatch, actual departure epoch/ownership, runtime source hooks and event accounting; external stationary sources retain the current upstream sampler. |
| `Exosphere_FiniteInventory.cpp` / `FiniteSurfaceRelease.h` / library makefile | Existing user-approved competing first-order depletion, weighted packets and finite remainder; required to avoid overdraw for residence times shorter than a particle step. |
| `Exosphere_Parallel.cpp` | Explicit prescribed/flux-balanced reservoir policy and signed, once-per-global-event external-supply diagnostics; prevents counting replicated MPI state several times. |
| `Exosphere_SourceProcesses.cpp` | Optional input-selected MIV rate/temperature with validated SI units; null hooks retain generated application constants. |
| `Exosphere_Chemistry.cpp` | Read-only actual product/loss notifications; no generic rate-law replacement. Moon's corrected Na dispatcher remains in the Moon application. |
| `src/pic/pic.cpp` / `.h` | Sampling-window reset independent of file-output flags; optional null-by-default event observer and exit classifier; finite departure-epoch scope is inert for legacy NaN calls. |
| `src/pic/pic_mover.cpp` | Optional moving-terrain route through the selected force mover, shortened collision-midpoint velocity and actual event notifications. Null geometry/observer retain normal upstream behavior. |

The current global runner's case-prefix protocol is preserved alongside Moon's
active-runner/log progress. Native Moon builds use bounded `make -j16` with
explicit `OPT3=-O3 -g`, after removal of the exact generated `build/` directory.
Root Makefile and site configuration are generated evidence, not source
changes for a commit. Detailed build/input/executable identities are retained.

Normal Moon builders and the registered Moon regression adapter now explicitly
select `-amps-test=off`. Config.pl intentionally preserves TESTMODE when only
the application name changes; after the canonical nightly target this otherwise
leaves a 100-step default cap plus nightly debug/sorting/load-balancing policy.
The initial delivered restore compiled successfully but inherited that mode;
its evidence is retained and superseded for normal-mode claims. The production
mode supporting campaign and focused registered Moon/harness retry recheck this
configuration correction. The canonical `test_Moon` target retains its explicit
`-amps-test=on`. Core Config.pl behavior and physical coefficients are unchanged.

## Regression acceptance

Run unmodified upstream and the integrated candidate with the same explicit
machine profile and existing verification fixtures. Compare exact suite/case
identities; every baseline PASS must remain PASS and no passing case may be
omitted. Existing upstream FAIL/ERROR/SKIP outcomes remain visible and cannot
be promoted by source preflight. The complete core registry has 1,085 cases at
the pinned revision; application outcomes have separate identities.

Existing missing native input assets and their recovery are described in
`../src/test/NATIVE_FIXTURE_PROVENANCE.md`. They are source fixtures, not model-
generated observations. Both baseline/candidate receive identical preserved
bytes. An initial incomplete-checkout attempt is retained separately.

Moon requalification uses the actual U-series and I01-I10 builds, followed by
one-/eight-rank sampling/event checks and the eight-rank M2 supporting campaign.
M1/I13 and formal M2 gates are unchanged. Upstream changes random-number
consumption and corrects source-angle sampling; a resulting nightly mismatch
must be investigated, with a separately reviewed candidate if needed. The
existing approved reference is not overwritten automatically.

## Current execution

The clean normal-mode Moon/SPICE build and all ten linked M0 checks now pass.
I04 initially had an evaluated numerical FAIL against its older reference:
four sampler values changed with the current-core source samplers. Two
four-rank runs of one executable and the delivered-checkout run produce
identical output. After the complete regression audit the user explicitly
approved the exact replacement SHA-256
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`
on November 12. It is installed with pinned reviewed provenance; the previous
reference, provenance and acceptance JSON are archived byte-for-byte. The
strict checker PASSes with a zero-byte diff and the completion message.
No tolerance or physical coefficient changed. No automatic replacement or
I13 freeze was performed.

Actual Na survival/event-ledger and sampling-window checks pass on one/eight
ranks. The eight-rank M2 supporting campaign has six local PASS, seven formal
I13-gated SKIPPED, six supporting PASS and one D02-readiness supporting
SKIPPED. The maximum Na/Ar relative inventory imbalance is about
`1.07e-15`. Formal M1/M2 qualification remains open. The complete upstream/
candidate registry comparison now PASSes preservation: all 1,193 cases are
accounted for, and every one of 919 baseline PASS results remains PASS. Both
sides have exactly 919 PASS, 272 SKIPPED, one FAIL and one ERROR, with zero
status changes or missing cases. The pre-existing nonpasses are
`T12-COVERAGE-STRICT-001` (open strict coverage) and
`CSWC-FULL-CSWC0628-MASS-TIME` (unchanged routine 180-second ceiling). This is
regression preservation, not an all-green registry or scientific validation.

The final application comparison exposed an additional infrastructure error
in Moon's existing I11-initialization reuse optimization: it queried approved
qualification grids before running CUSTOM regression I12. The fix restricts
receipt reuse to qualification mode. The real eight-rank MOON-I12 rerun then
evaluates unchanged production wiring, retains receipt/stdout/stderr/input/
manifest, and returns its original unqualified-CUSTOM SKIPPED. The first ERROR
attempt remains retained; no threshold, profile membership or physical input
was changed. The twelve campaign-infrastructure protocol checks pass. The
explicit normal-mode Moon/harness retry has 27 PASS and 22 SKIPPED, zero
FAIL/ERROR; the final whole-registry audit includes that exact retry.

The separate core and application checkouts had byte-identical production
source across all 653 compared C/C++/Fortran/header/definition/library-makefile
paths. Later changes are documented Python control-plane fixes and README/
status updates, with focused checks. The delivered-checkout nightly and
canonical Moon/SPICE rebuild are recorded separately from the frozen complete
registry campaign. Normal-mode Na/Ar, budget and sampling rechecks PASS after
the configuration correction. The normal Na-only Moon/SPICE executable is
restored with `TESTMODE=off`, `-O3 -g` and one OpenMP thread per MPI rank.

Durable evidence is in
`test_output/srcMoon/core-integration-20261112T224309Z/`. Its completion report
and `results/regression-comparison-final/{result.json,comparison.csv}` retain
case membership, counts, source/data/build hashes, exact commands, known
nonpasses and all initial attempts/retries. Original absolute `/tmp` paths in
raw receipts are preserved; the evidence-copy README maps them to the durable
tree. Backups and reference-review history are retained there as well.

## Reproducing the regression comparison

Use separate checkouts for upstream and candidate. Do not configure, build or
run another native owner in a checkout already owned by a campaign. Each local
checkout needs its own installed SHARE library; complete `make -j4 share_lib`
before starting the native suites. The recorded `machine.json` preserves the
maintained Passion tool/capability profile and relocates only work/output/data
roots. It selects GCC 9.4, OpenMPI 4.0.3 and the existing SPICE installation.
No test threshold or native owner timeout is changed.

For the two generated Enceladus receipts, prepare the path metadata in each
application checkout before its campaign:

```sh
python3 -B srcMoon/test/prepare_core_regression_inputs.py \
  --repo-root "$PWD" --output-dir /tmp/AMPS-NEW-input-metadata
```

Every line other than the three named absolute-path comments must remain
byte-identical. The existing exporter and test owner remain unchanged. The
same helper/fixtures prepare upstream and candidate; retained initial setup
errors are not erased by the later valid reruns.

From each appropriately prepared checkout, with its own fresh output path:

```sh
python3 -B src/test/run_tests.py --config /path/to/machine.json \
  --match 'CORE-*' --streams 1 --jobs 8 --build-jobs 16 \
  --output-dir /tmp/AMPS-NEW-core-results
python3 -B src/test/run_tests.py --config /path/to/machine.json \
  --match 'APPLICATION-*' --streams 2 --jobs 2 --build-jobs 16 \
  --output-dir /tmp/AMPS-NEW-application-results
```

The complete pinned registry selects 1,085 core and 108 application cases.
An optional application preflight SKIP is not a linked Moon PASS; separate
Moon M0/supporting runs provide the production evidence above. Native
Enceladus E3 retains its registered one-rank, 10,000-step contract. Moon
supporting runs use eight MPI ranks and one OpenMP thread per rank.

Audit the exact aggregate reports after execution (add each justified retry
as another argument in explicit attempt order):

```sh
python3 -B srcMoon/test/compare_core_regressions.py \
  --baseline /path/to/upstream-core/amps-core-tests.json \
  --baseline /path/to/upstream-applications/amps-core-tests.json \
  --candidate /path/to/candidate-core/amps-core-tests.json \
  --candidate /path/to/candidate-applications/amps-core-tests.json \
  --output-dir /tmp/AMPS-NEW-regression-audit
```

Each aggregate needs its original neighboring `selection.json`. The audit
keeps every actual earlier baseline PASS protected even when a later retry
fails. Missing rows become ERROR/SKIP according to their infrastructure
evidence, never inherited PASS. It records exact report hashes, all status
changes, and a per-case CSV with runner/log locations. A preservation PASS
means no passing regression was lost; existing upstream FAIL/ERROR remain
visible. Three protocol checks run automatically before the audit. They test
report handling, not scientific observations.
