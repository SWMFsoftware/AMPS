# Global test-runner governance and rejection controls

Updated 2026-11-10. `CORE-HARNESS-001` owns 13 cases. It tests the real runner
and its read-only validation authorities, using isolated synthetic fixtures
where deliberate bad input/process behavior is needed. A harness PASS proves
that control-plane contract; it does not qualify production physics or upgrade
partial source coverage.

## Run the owner or the new coverage case

From the AMPS root:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test CORE-HARNESS-001 --output-dir test_output/harness-check

python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --case HARNESS-COVERAGE-GATES-001 --output-dir test_output/coverage-gate-check
```

Choose a fresh output directory. Another host must use a real profile based on
`test/qualification-machine/template.json`; `--config` is required. Nested
commands receive explicit portable Python/profile argv and local work paths.
Tests do not embed the author's compiler, home, output or dataset paths.
The harness does not rebuild native AMPS. Native-configuration scenarios use
private synthetic tools to exercise the maintained adapter control flow.

## Complete case inventory

| Case | What / how / why |
| --- | --- |
| `HARNESS-CONFIG-001` | Layered defaults, machine files, environment and CLI precedence, path anchoring and required configuration; prevents a worker silently using another site's paths/options. |
| `HARNESS-NATIVE-CONFIG-ISOLATION-001` | Actual native adapters with private valid and deliberately contradictory saved/effective no-SPICE configurations and build receipts; prevents stale configuration or link flags from reaching MPI as an apparent physics failure. Implementation: `native_configuration_isolation.py`. |
| `HARNESS-MANIFEST-001` | Manifest/component inventory consistency and rejected malformed/missing registrations; protects exact case ownership. |
| `HARNESS-SELECTION-001` | Exact suite/case/campaign/alias selection and rejected invalid requests; protects what the report claims was selected. |
| `HARNESS-PROTOCOL-001` | Valid and hostile result records, counts/status/exit consistency and complete case membership; prevents hidden missing results or a smaller green report. |
| `HARNESS-SCHEDULER-001` | Concurrent scheduling plus exclusive barriers; protects native checkout-mutating operations from overlapping eligible jobs. |
| `HARNESS-REPORT-001` | Final individual rows and aggregate/stream count reconciliation; every required case retains a visible state. |
| `HARNESS-SCHEMA-001` | Schema version and campaign/capability definition rejection; prevents unsupported settings being silently accepted. |
| `HARNESS-SUBPROCESS-001` | Hostile child receipts and subprocess failures; distinguishes execution/control errors from successful assertions. |
| `HARNESS-TIMEOUT-001` | Timeout/interruption process-group cleanup, including descendants; prevents abandoned test processes continuing after a reported timeout. |
| `HARNESS-CONCURRENCY-001` | Bounded job counts, completion/fail-fast scheduling and retained results for already-started work; prevents the CLI concurrency request from being only a label. |
| `HARNESS-TRACEABILITY-001` | Expanded global case ownership and source/configuration traceability; prevents stale or unregistered identities being promoted to coverage. |
| `HARNESS-COVERAGE-GATES-001` | Actual standalone/strict CLIs, contradictory declarations and an injected actual-CLI governance defect; protects the strict completeness gate. Implementation: `coverage_gates.py`. |

The main owner is `run_test.py`. The two named helper files isolate substantial
configuration and coverage scenarios; the global manifest remains the single
registration authority.

## New strict-coverage contract in detail

`validate_traceability.audit_repository()` is now the shared read-only
source/configuration/catalog ownership audit used by both the standalone
validator and the global runner's strict option. The global strict path must
not rely solely on component labels while contradictory source ledgers remain.

The new case executes:

1. The actual `validate_traceability.py --output ...` CLI and retains expanded
   source ownership, nonzero translation-unit/case counts and its exit 0.
2. Private inventory controls: a synthetically complete fixture is accepted
   only as a fixture; partial coverage remains incomplete; missing/unknown
   owners and a newly discovered unregistered component reject.
3. The actual owner-validation function: `tested` with required gaps, no direct
   cases or stale case IDs rejects. Synthetic declarations never modify the
   real inventory or supply scientific evidence.
4. The real `run_tests.py --coverage --require-complete-coverage` CLI with an
   explicit generated profile. A valid incomplete real inventory must exit 1;
   only an actually complete, consistent inventory can exit 0.
5. A private child-process fault injection at the real validator interface,
   then execution of the unmodified real strict CLI. It must return exit 2
   containing `INJECTED-COVERAGE-LEDGER-CONTRADICTION-001`. This catches bypassing
   the cross-ledger audit even when ordinary incomplete coverage would still
   return 1. No production file is patched for this control.

There are ten retained scenario descriptions under this one registered case.
The distinction is essential:

| Condition | Strict CLI exit | Harness interpretation | T12 interpretation |
| --- | ---: | --- | --- |
| Consistent real inventory with incomplete components | 1 | PASS if this truthful result is preserved | `T12-COVERAGE-STRICT-001` FAIL |
| Contradictory source/configuration/catalog ownership or broken invocation | 2 | PASS only for the deliberately injected/rejection control; otherwise harness failure | ERROR, not scientific FAIL |
| Actually complete, consistent coverage | 0 | PASS only if the real inventory supports it | Coverage gate can PASS |

A PASS of the TEST OF THE GATE does not turn the incomplete gate into PASS.
The current real inventory has 2 complete and 26 incomplete components. T12
continues to report that FAIL. No coverage label was upgraded by these changes.

## Evidence, failures and current result

The global suite row records its log, artifact and work paths. Coverage-case
artifacts include:

```text
standalone-traceability.txt
standalone-traceability.json
strict-coverage-cli.txt
strict-coverage-injected-defect.txt
coverage-gate-controls.json
```

The JSON retains scenario descriptions/count, actual CLI exits and the real
incomplete-component list. The generated private profile/bootstrap remain in
the work tree when the machine's keep-work policy retains it. Required malformed
fixtures are not omitted because a platform/tool is inconvenient.

The complete 2026-11-10 core campaign passes all 13 harness cases. Focused
`harness-final-controls-corrected/amps-core-tests.json` also passes 13/13 with
the injected-defect control. Evidence root:
`test_output/remaining-failures-repair-20261110/`. First schema/fixture formatting
failures and their correction remain retained separately; no real inventory or
strict exit was changed to hide them.

## Contributor / agent rule and limits

New cases must update the manifest/registry, source/configuration ownership,
this README and generated `src/test/TEST_CATALOG.md` in the same change. Run
catalog generation/check and standalone traceability, then execute the real
owner and inspect every case state. A fault fixture is useful only if the
wrong runner behavior makes the control fail; a test-local replacement must
not be mistaken for production execution.

These controls do not close all cross-platform, cross-stream, cancellation,
sharding, filesystem-isolation, legacy-oracle or resource-stress gaps. The
complete core result and remaining physics/coverage work are in
[`src/test/STATUS.md`](../../../STATUS.md) and
[the global procedure](../../../README.md).
