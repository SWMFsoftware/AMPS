# Repair plan: native exosphere and dust configuration isolation

Prepared 2026-11-10. Implementation update: the runner repairs, shared mode
verification, new harness case and documentation are implemented. The complete
harness and real starting-state/reverse-order native matrix pass. Full core
baseline revalidation PASS: all 785 original passes preserved, exactly the two
requested failures repaired, and seven new cases PASS. Final counts are
1,035 total / 792 PASS / 2 FAIL / 240 SKIP / 1 ERROR, with the original
nonpasses unchanged. Exact receipts are linked from STATUS.md and retained
under `test_output/native-config-isolation-20261110/` as local evidence.
The original execution plan below is retained as the acceptance checklist.

## Goal and acceptance boundary

Restore these existing cases to PASS regardless of whether the checkout was
previously configured with a real SPICE installation or with `nospice`:

| Owner | Existing individual case |
| --- | --- |
| `CORE-EXOSPHERE-SPHERE-NATIVE-001` | `CORE-EXOSPHERE-SPHERE-INTEGRATION-001` |
| `CORE-DUST-NATIVE-001` | `CORE-DUST-CHARGING-DYNAMICS-001` |

Keep their case identities, physical inputs, numerical oracles, tolerances,
one-rank execution, exclusive checkout ownership, and source/reference data.
This repair concerns test build configuration. It does not implement geometry
acceleration or close the independent 3-D Hartley2 runtime/physical gates.

## Evidence and cause

Failing run:
`/data/vtenishe/AMPS-test-results/run-20261110T092615Z-924820/`.

- Exosphere builds but aborts in `utc2et_c` with `SPICE(NOLEAPSECONDS)`.
  Evidence: `artifacts/CORE-EXOSPHERE-SPHERE-NATIVE-001/exosphere-native-output.txt`.
- Dust fails to link with undefined `Exosphere::IAU_FRAME`.
  Evidence: `artifacts/CORE-DUST-NATIVE-001/dust-native-build.txt`.
- The earlier passing exosphere build contains `_NO_SPICE_CALLS_` and does not
  link CSPICE. The failing build omits that compiler define and links CSPICE.
- Both adapters call Config.pl with only an application argument. Config.pl
  preserves the global SPICE installation unless explicitly changed through
  `-spice-path=...`. The component's input-block `SPICE=off` is insufficient to
  establish the global build mode.
- The established top-level Makefile supplies `_NO_SPICE_CALLS_` when its
  effective `SPICE` value is `nospice` or empty. Makefile.local is included
  after the defaults and can retain an application's real installation path.

The original baseline contained two FAIL cases and one ERROR case, while only
one suite was classified FAIL. A failed case also belongs to the ERROR suite.
The new run adds exactly two PASS-to-FAIL changes; current baseline retention
is 783/785. The six newly registered comet/input cases still PASS.

## Step 1: pin the two fixture configurations

Modify only the test adapters initially:

- `src/test/tests/models/exosphere_native/run_test.py`
- `src/test/tests/models/dust_native/run_test.py`

Each must select its application and explicitly pass `-spice-path=nospice`
through Config.pl before compiling. This can be one structured argv list:

```text
./Config.pl -application=test/core-exosphere-sphere -spice-path=nospice
./Config.pl -application=test/model-dust -spice-path=nospice
```

Retain the exact existing native preflight for every rebuild: confirm the
AMPS root and idle native state; acquire the existing exclusive lock; remove
only `build/`; regenerate the selected application; compile through top-level
make with `MAKEFLAGS=-j16`; retain and run the executable on one MPI rank.
Preserve the global stream gate and hold exclusive ownership through the
configuration, build, verification, and runtime phases.

Do not edit `.amps.conf` directly, change Config.pl's global default, disable
SPICE globally in the orchestrator, add dummy frame definitions, or load new
kernels into these fixtures. SPICE-dependent application owners must remain
able to select their own installation normally. Keep Makefile.conf and site
libraries intact; continue using machine-profile tools and paths.

## Step 2: verify effective configuration and preserve evidence

Add a small shared test helper, proposed location
`src/test/lib/native_configuration.py`, for read-only mode inspection and
consistent diagnostics. Keep application selection explicit in each owner;
this helper is not a rewrite of the native build framework.

Before compiling, verify:

1. The selected application/input is the intended fixture.
2. The effective global `SPICE` mode is `nospice`.
3. The actual make compiler options define `_NO_SPICE_CALLS_`.

Inspect both the saved selection and the effective make settings. Checking
only for the text `SPICE=nospice` in one file is insufficient when multiple
make inputs or command-line overrides can disagree. A temporary probe makefile
can add a dedicated, dependency-free inspection target and print `SPICE` and
`SEARCH_C` after the real Makefile is read. Invoke only that target using the
same configured make tool/environment; do not use `make -n amps`, whose
recursive make recipes can still execute. Keep the probe in work/artifact
storage, not in the production Makefile or Git control files.

Record the observed values and pass/fail results in a retained configuration
receipt. Update build-identity JSON to contain the exact new configure argv,
expected/observed mode, executable/deck hashes, and relevant configuration
hashes. Retain configure output even if configuration or verification fails.
After compilation retain/check the actual build commands as corroborating
evidence that the intended define was used and CSPICE was not selected.

A contradictory or unverified mode must produce a clear ERROR for the owned
case before native execution, with the relevant log/receipt location. It must
not become PASS, SKIP, or an unrelated numerical failure. Keep all existing
native source-ledger and charging-oracle acceptance criteria unchanged.

## Step 3: durable control-plane regression tests

Propose one additional case under the existing `CORE-HARNESS-001` owner:
`HARNESS-NATIVE-CONFIG-ISOLATION-001`.

Exercise the adapters' control flow with isolated external-tool fixtures:

- Starting mode enabled versus disabled: both fixtures explicitly establish
  the required mode before the builder is invoked.
- A Config.pl stand-in returns success but ignores the requested mode: the
  adapter rejects the mismatch and never starts native execution.
- Saved selection and effective make mode disagree: reject the inconsistency.
- Configuration command failure: retain its output and produce a nonpassing
  result with the correct case identity.
- Build-identity/configuration receipts describe the actual invocation and
  observed mode, including the newly required argument.

These are synthetic harness checks of configuration and reporting, not physics
qualification. Mock native receipts must remain confined to the harness's
private fixtures and must never count as executions of the real physics cases.
Actual native runs in Step 4 provide the independent end-to-end proof.

Register the new case in the manifest/harness and appropriate configuration
contract/traceability records. Update its what/how/why/where description and
regenerate `TEST_CATALOG.md`, following the agent catalog rule. The proposed
case is not registered by this plan. Existing physics case IDs remain intact.

## Step 4: native validation from both starting states

Use the machine profile to obtain a real SPICE path when available; do not
hard-code Passion paths in the adapters. Use Config.pl to seed starting state,
with the checkout idle. No Hartley2 compilation or long initialization is
needed merely to establish a SPICE-enabled starting configuration.

Reset the starting state independently before every row:

| Initial global selection | Fixture | Required result |
| --- | --- | --- |
| SPICE-enabled application | Exosphere | PASS, verified no-SPICE build |
| SPICE-enabled application | Dust charging | PASS, verified no-SPICE build |
| `nospice` | Exosphere | PASS, same declared profile and numerical checks |
| `nospice` | Dust charging | PASS, same declared profile and numerical checks |

The second row is essential: running dust only after the repaired exosphere
owner could conceal its own failure to reset configuration. Also run the pair
in reverse suite order from an enabled start. Each native row uses a clean
build, one MPI process, private artifacts, and an unused output directory.
Retain executable hashes independently; require consistent generated modes
and numerical invariants rather than assuming separate builds are byte-identical.

Verify the unchanged physical receipts:

- Exosphere: one exact PASS marker, positive counts for all five compiled
  source modes, AMR/particle ownership and represented-rate/inventory closure.
- Dust: one exact PASS marker, all four radii, finite charging histories,
  existing published-equation current/root checks and existing convergence
  and `q/a` tolerances.

Finally select a legitimate SPICE-enabled application through Config.pl and
inspect its mode to confirm that the repair has not changed normal SPICE
selection. This is a configuration boundary check, not a completed Hartley2
runtime qualification.

## Step 5: audit adjacent owners and run the full regression

Audit other application-only native selectors, particularly
`tests/pic/ecsim_operators/run_test.py` and
`tests/models/dust_injection_native/run_test.py`. Pin the same no-SPICE mode
only where the fixture's declared contract requires it. Native core and mesh
owners already explicitly disable SPICE and provide useful consistency
examples. Preserve the profile-driven SPICE requirements of application tests.
Document other inherited settings found by the audit; do not silently change
unrelated compile profiles without an established fixture requirement.

After focused runs pass, seed a SPICE-enabled start and run the complete core
selection using the same Passion profile and concurrency as the reported run.
Compare individual case identities with both the frozen pre-integration
baseline and the failing user run. Require:

- Both newly failing cases become PASS.
- All original 785 baseline PASS identities remain PASS; no new missing/SKIP
  cases are accepted as preservation.
- All six added comet/input cases and all new harness checks PASS.
- `CSWC-FULL-CMBGU05`, `T12-COVERAGE-STRICT-001`, and
  `CSWC-FULL-CSWC0628-MASS-TIME` remain separately reported with their existing
  FAIL/FAIL/ERROR identities; they are outside this repair.
- Every other baseline status remains accounted for. Any unexpected change
  requires investigation before acceptance.

For the original 1,034 identities the restored target is 791 PASS / 2 FAIL /
240 SKIP / 1 ERROR. With exactly the one proposed additional harness case, the
core target becomes 1,035 total / 792 PASS / 2 FAIL / 240 SKIP / 1 ERROR. The
suite inventory remains 50 (43 PASS / 1 FAIL / 5 SKIP / 1 ERROR). Reconcile all
generated inventory counts when registering the case; current counts do not
change merely because this plan proposes a test. The expected overall runner
exit remains nonzero because the unrelated existing failures/error persist.

Commands after implementation, from the AMPS root (choose fresh output paths):

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-HARNESS-001 --output-dir test_output/native-config-isolation/harness

python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-EXOSPHERE-SPHERE-NATIVE-001 --test CORE-DUST-NATIVE-001 --streams 1 --jobs 1 --output-dir test_output/native-config-isolation/pair-normal

python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-EXOSPHERE-SPHERE-NATIVE-001 --test CORE-DUST-NATIVE-001 --order reverse --streams 1 --jobs 1 --output-dir test_output/native-config-isolation/pair-reverse

python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --match 'CORE-*' --streams 1 --jobs 8 --output-dir test_output/native-config-isolation/core-final

python3 -B src/test/generate_test_catalog.py --check
python3 -B src/test/validate_traceability.py
```

The commands alone do not seed the four starting-state matrix rows; retain
those separate Config.pl setup commands and their observed settings with each
row's evidence. On other machines substitute the required machine profile.

## Step 6: documentation and completion record

Update the two fixture READMEs and the native-testing documentation with the
explicit build-mode contract, manual Config.pl commands, evidence fields, and
failure diagnostics. Keep detailed comments explaining global build selection
versus component input settings. Update `src/test/STATUS.md` and
`COMET_CORE_INTEGRATION_STATUS.md` with exact commands, fresh report paths,
PASS/FAIL/SKIP/ERROR counts and the baseline comparison.

Generated `test_output/`, build products and current-run logs remain local
evidence, excluded from the source import. Preserve all existing reports,
user-owned changes, site files and physical inputs. Git remains read-only;
leave every source/documentation change unstaged and uncommitted.

Completion requires the real native matrix and full baseline comparison,
not just adding a configure flag or passing the synthetic harness case.
