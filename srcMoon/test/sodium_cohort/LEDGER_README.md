# Independent event-ledger verifier

Run the maintained controls from the AMPS checkout:

```bash
python3 -B srcMoon/test/test_m3_sodium_ledger.py --output-dir /absolute/new/output
```

This bounded Python component imports no production mass/kinematic/chemistry
formula. It reports PASS/FAIL/ERROR/SKIPPED with exits0/1/2/77, retains unittest
logs/source hashes and generates CSV/PNG/EPS for a labelled arithmetic toy
ledger. These figures test bookkeeping software; they contain no observations
and do not qualify physical sodium kinetics. Missing retained real source
packets are separately skipped. The current actual partial source profile
refuses unavailable Na+ before event scoring. No source value, threshold,
branch population, product state or tolerance is supplied to fill that gap.

`verify_ledger.verify(paths,profile_json,mapping,bounds)` consumes per-rank
`event_metrics.rankR.hex.csv` with exact columns:

```text
original_parent_id,key,value
```

`value` must be a finite native hexadecimal binary floating representation.
The parser converts its significand/exponent directly to `fractions.Fraction`,
preserving long-double bits and full exponent range. Decimal source CSV values
are also exact rationals. There is no intermediate Python float in ledger,
conservation, Q, geometry-vector or native-readback comparisons. The only
binary64 conversion checks the documented PIC storage of the exact SI charge.

Mapping contains `parent`, `ion`, `electron` source state IDs,
`profile_sha256`, `event_frame`, `ion_compiled_species` and
`ion_catalog_species`. These are explicit declarations from the frozen native
input/compile/profile receipts; they are never inferred from species names.
A source-bound experimental run also requires `photon_source`:

```json
{
  "reaction_id": "PHIDRATES:Na:00:na_e",
  "solar_activity": "0",
  "rates": {"path":".../photo_rates_1au.csv","sha256":"..."},
  "thresholds": {"path":".../thresholds.csv","sha256":"..."}
}
```

Only a unique native quiet/active endpoint at1AU, optical_depth0 is selected.
The native endpoint rate, mean excess energy and threshold are independently
read/authenticated and compared to the serialized source receipt; no source
mean defines a complete daughter PDF or repairs a state/threshold mismatch.

Bounds use schema `srcMoon-M3-event-ledger-bounds-v1`, explicit `event_frame`,
`scope` (`source-bound-experimental` or `analytical-software-fixture`) and
`mass_kg`, `momentum_kg_m_s`, `energy_j`, `number`, `rotation` objects, each with
`absolute`/`relative` decimal strings. `rate_s` is additionally required for a
source endpoint comparison. The names `mass`, `momentum`, `energy` are accepted
aliases for existing frozen planner bounds. There are no numerical defaults.
Frame-velocity residuals are mass-weighted before applying momentum bounds,
so m/s and kg*m/s tolerances are never confused. Number bounds govern statistical
weight round-trip arithmetic and remain distinct from survival confidence or
numerical-bias allocations. Caller freezes these bounds and authenticates the
actual native producer/input/build before invocation; a positive verifier does
not establish that independent admission. Existing0.2/0.3 percent bias and
formal ensemble limits stay in the qualification scorer.

The verifier independently reconstructs source masses, formation Q, integer
nuclei/charge, represented stoichiometry and both per-event and represented
momentum/kinetic+state energy. It assembles the reservoir from the outgoing
free electron minus incoming photon, then applies the actual proper rotation
and velocity offset; it never defines the reservoir as the negative measured
tracked residual. Photon masslessness and parent-rest Doppler energy are
checked through exact squared interval inequalities without square-root
rounding. Actual committed initial PIC child velocity/position/weight/internal
energy/species/state/epoch/original interval must equal its immutable plan
EXACTLY. Actual compiled PIC masses and SI charges are checked independently
against the authenticated source profile. Native Boolean verified flags are
ignored. Empty/missing capability is SKIPPED; malformed/duplicate/changed input
is ERROR; an evaluated mismatch is FAIL. A later missing input cannot hide an
already evaluated FAIL. See [LEDGER_CONTRACT.md](LEDGER_CONTRACT.md).

The main native runner retains decimal21 metrics for inspection and exact hex
metrics for this verifier. Its source/manifest must bind the actual getter
readback immediately after successful transaction staging and before child
continuation, not a later evolved child. This verifier cannot turn an arbitrary
CSV into actual runtime evidence: native executable/source/configuration,
predecessor, numerical-convergence and frozen-protocol bindings remain required.

The return value verifies the events supplied, not completeness of a cohort
selection. The owning conditional executor must reconcile the set of original
parent IDs in all native hex metric files against every REACTED entry in the
actual `parent_records.csv`, bind every rank's artifacts/native completion, and
reject omitted/duplicate records before using the result. It must retain the
actual frozen input/protocol hashes and source-state/operator preflight. This
component never treats a scalar actual/reference JSON or subset of events as a
complete independent convergence proof.
