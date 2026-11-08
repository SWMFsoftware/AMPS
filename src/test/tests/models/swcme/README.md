# SWCME model-local registry adapter

This adapter imports the SWCME-owned test surface into the portable global
harness while keeping three different qualification scopes separate.

## Imported scopes

1. The bounded Phase-R1 runner in
   `src/models/swcme/test/run_tests.py` owns four executable cases. The adapter
   copies the complete `src/models/swcme` tree into its per-run work directory,
   omits prior build products, and executes that copied runner with `--all`.
   The local JSON report is checked for its exact schema, ordered membership,
   totals, unique IDs, statuses, and exit-code agreement before its results are
   mapped to global IDs.
2. `src/models/swcme/test/core/test_swcme.cpp` retains 127 extended cases, but
   the Makefile currently names six unavailable translation units. The adapter
   validates the complete C++ registry and the exact missing-source set, then
   reports every extended case as `SKIP`. It never treats R1 as evidence for
   this wider registry.
3. `src/models/swcme/test/validation/campaign.json` retains VP01--VP16 in
   dependency order. Each entry and its case-owned `case.json` are validated.
   The authoritative `validation/run_validation.py` is absent, so root cases
   are reported as typed `SKIP` and dependency-bearing cases as typed
   `BLOCKED`, represented through the global protocol's `SKIP` status and an
   explicit local-status message. No observational qualification is claimed.
   VP01's `case.json` also names an absent `run_vp01.py`; that additional
   limitation is retained in the audit and VP01 result message.

The complete global inventory is 147 cases: 4 R1, 127 extended, and 16
observational. Registry validation always precedes selected-case filtering.
If a missing extended source or the observational runner is restored, the
adapter intentionally returns `ERROR`; its old non-executable classification
must then be replaced with an executable adapter rather than silently retained.

## Result and evidence contract

The adapter writes the normal schema-1 global result. Evidence includes:

- `swcme-registry-audit.json`, containing source hashes, exact local IDs,
  counts, and current blockers;
- `swcme-r1-transcript.txt` when any R1 case is selected;
- the authoritative R1 JSON and JUnit reports when R1 reaches reporting.

Suite status precedence is `ERROR`, `FAIL`, `SKIP`, then `PASS`. Therefore a
full 147-case selection is necessarily `SKIP` at the current source boundary,
even when all four executable R1 cases pass. Selecting only passing R1 cases
can return `PASS`, but still performs a complete registry-drift audit first.

## Global IDs

- R1: `SWCME-R1CFG01`, `SWCME-R1D01`, `SWCME-R3D01`, `SWCME-R13D01`.
- Extended: `SWCME-EXT-<local-id>`. The legacy `1D_AMBIENT_01` spelling is
  normalized to `SWCME-EXT-1D-AMBIENT-01` to satisfy the global ID grammar.
- Observational: `SWCME-VP01` through `SWCME-VP16`.

## Direct synthetic-context test

The adapter uses the same context contract as every `src/test` runner. A
standalone check should create a schema-1 context whose `expected_cases` is the
complete ordered `ALL_GLOBAL_CASES` tuple and whose `selected_cases` is the
desired subset, then invoke:

```sh
python3 src/test/tests/models/swcme/run_test.py --context /tmp/context.json
```

Required configured tools are `python` and, when R1 is selected, the local
runner's normal `make`/C++17 toolchain. The adapter does not perform a native
AMPS application configuration or checkout-wide build.
