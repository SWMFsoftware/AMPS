# Audit a completed sodium commissioning run

This is a thin independent supporting audit of actual artifacts, not an AMPS
launcher. Read [the contract](COMMISSIONING_AUDIT_CONTRACT.md) first. It neither
builds nor advances a particle, modifies source data, selects a budget, or
creates a qualification plan. The P02 mechanical-cap audit is separate.

Use the same original-owner binding as `verify_p02.py`:

```json
{
  "schema": "srcMoon-M3-P02-owner-binding-v1",
  "profile": {"path": "/absolute/profile.json", "sha256": "ACTUAL_SHA256"},
  "budget": {"path": "/absolute/approved-budget.json", "sha256": "ACTUAL_SHA256"},
  "owners": [{
    "owner_manifest": {"path": "/completed/run/run_manifest.json", "sha256": "ORIGINAL_SHA256"},
    "owner_result": {"path": "/completed/run/result.json", "sha256": "ORIGINAL_SHA256"},
    "mpi_ranks": 8
  }]
}
```

The placeholders above are intentionally not valid descriptors. Pin the real
completed wrapper receipts immediately after the original native owner returns.
The manifest's original artifact map authenticates the native input receipt,
parent CSV and every rank's exact hexadecimal event CSV. A newly rehashed CSV
alone cannot replace original owner evidence. A source-only/preflight run, a
failed owner, or an arbitrary Boolean cannot qualify as a completed physical
owner. Use previously absent output directories:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/audit_commissioning.py \
  --binding /absolute/original-owner-binding.json \
  --output-dir /absolute/new-commissioning-audit
```

Expect `result.json`, `run_manifest.json`, `audit.log`, `comparison.csv`,
`comparison.png`, `comparison.eps`, `survival_diagnostics.json`, and per-owner
`event_time_links-XX.json`, `inverse_square-XX.json`, and `ledger-XX/` artifacts.
Ledger directories retain the maintained independent rational verifier's own
CSV/PNG/EPS, result and manifest. Failures and partial evidence are preserved.
Matplotlib is required for figures; missing plotting dependencies produce an
explicit SKIPPED component while retaining already evaluated criteria.

The main status/exit code is PASS/0, FAIL/1, ERROR/2, or SKIPPED/77 for the
declared supporting checks. It always separately reports:

- `formal_M3_status = SKIPPED / NOT QUALIFIED`;
- `refinement_status = NOT_EVALUATED`;
- `statistical_qualification_status = NOT_EVALUATED`;
- `P02_cap_audit_status = SEPARATE_OWNER_NOT_EVALUATED_HERE`.

The independent clock reconstructs each original parent's semantic HMAC source
uniform and analytical `-ln(U)/lambda`, using the actual source normalization,
visibility and geometry from the authenticated receipt and the retained native
D04 endpoint. Source tables/QA/provenance are authenticated by the existing
source parser. The existing **0.002** normalized first-hit timing allocation is
checked at this one level; this does not prove refinement independence. The
fixed **0.001** I21 inverse-square criterion is checked independently for actual
accepted rates. No event means that rate comparison is SKIPPED, not a fabricated
rate PASS. The complete parent range must agree exactly with the per-rank event
set, and loss through a surface/domain is never silently censored.

Survival diagnostics show the actual survivor/removal counts and independent
expected survival. There is no lifetime/branch acceptance test, confidence claim,
fixed-N adequacy assertion or agreement-based selection. A 128-parent run can
demonstrate actual source/transaction wiring and event bookkeeping while
remaining insufficient for the frozen I21 ensemble/statistical campaign.

Run software controls without native data:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/test_audit_commissioning.py
```

These tests use analytical exponential/parser fixtures and missing/corrupt
administrative bindings only. They never fabricate a successful native owner.
