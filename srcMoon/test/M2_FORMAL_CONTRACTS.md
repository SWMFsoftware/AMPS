# Shared M2 formal-contract dispatch

The ordinary `.in` setting `m2_formal_contract_file` can identify a single
I15/I19/I20 document, or a reviewed `srcMoon-M2-formal-contract-set-v1` bundle.
This permits one preserved application receipt to dispatch different external
contracts for illumination, the complete native numerical campaign and polar
comparison. These meanings are not merged and their predecessors are not
waived. The direct document still undergoes its owning test's schema/scientific
checks. A set approval does not approve an unreviewed child or fabricate data.

The bundle contains `schema`, `status=APPROVED`, `review`, and `contracts`.
`review` requires literal `approved=true`, a nonempty `reviewer`, and an
`approved_utc` including its time zone. Each selected `contracts.I15`,
`contracts.I19` or `contracts.I20` record is exactly a `path` and a lowercase
64-digit `sha256`. The path resolves against the bundle directory. Every child
must exist and its actual bytes must match. Nested bundles and self-reference
are rejected; the resolver never searches fallback locations or guesses a
child from a filename.

Child schemas remain:

| Test | Direct schema / documentation |
|---|---|
| I15 | `srcMoon-I15-reference-contract-v1`; [I15 README](integration-tests/I15/README.md) |
| I19 | `srcMoon-I19-formal-execution-v1`; [I19 README](integration-tests/I19/README.md) |
| I20 | `srcMoon-I20-polar-contract-v1`; [I20 README](integration-tests/I20/README.md) |

The run manifest preserves both bundle and selected-child path/hash records.
There is no automatic freeze, reference installation or default READY bundle.
The investigator must supply actual approved files. Missing set/child/review
or test membership is SKIPPED; malformed/ambiguous JSON or changed hashes is
ERROR. A hash-valid bundle is dispatch provenance only, never a linked or
observational PASS. U/I IDs and scientific thresholds are unchanged.

The metadata-only refusal and dispatch tests run with:

```bash
python3 -B srcMoon/test/test_m2_contracts.py
```

Fixtures contain no observations or physical model values, remain temporary,
and are not qualification evidence. Tests cover relative/case/space-preserving
paths, both document hashes, missing test IDs, caller-specific unavailable
status, absent review, duplicate/non-finite JSON, bad hashes, changed children,
nesting/container corruption and explicit review time zones.
