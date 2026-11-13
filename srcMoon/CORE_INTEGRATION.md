# Integration with the current AMPS core

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
