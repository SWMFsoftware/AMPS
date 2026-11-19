# Integrated exposure and optical-clock bridge

Supporting coverage extends the existing U15/U16 reaction mechanics; no new
permanent ID is assigned. The contract is
[`M3_INDEPENDENT_SOFTWARE_CONTRACT.md`](../../M3_INDEPENDENT_SOFTWARE_CONTRACT.md).
The normal Moon library compiles `ReactionExposure.cpp`, `PiecewiseHazard.cpp`
and `ReactionExecution.cpp`. This runner links those exact production files.

```bash
python3 -B srcMoon/test/reaction_execution/exposure_bridge/run.py \
  --output-dir test_output/srcMoon/exposure-bridge-NEW
```

Rates are explicitly declared affine functions of relative physical seconds;
exposure is dimensionless. They are mathematical controls, not observations,
solar inputs or a permitted approximation of a lunar trajectory. The main
oracle is `H(t)=t+t²`, with closed-form first crossing and instantaneous
branch probabilities. Its integrated fractions deliberately select a different
channel, detecting an otherwise plausible process-selection error.

The runner also checks immutable preview/retry, compensated small debits,
partitioning, gaps (including future gaps beyond an early crossing), shadows,
consumed clocks, invalid draws and left-limit discontinuity semantics. A
crossing with zero instantaneous rate is explicitly refused; the service does
not move the event or invent a probability. Independent analytical CSV, PNG
and EPS diagnostics are preserved with compiler commands, source/executable
hashes, stdout/stderr, manifest and result. PASS is local numerical/mechanical
verification only. Missing compiler/plot dependency is SKIPPED; infrastructure
errors are ERROR; violated evaluated identities are FAIL. Exit codes are
0/1/2/77. No RNG is generated here; supplied open-uniform controls are explicit.

The trajectory-dependent controller now consumes the same full-component
bridge through an explicitly installed accepted-arc provider. Physical
product laws and complete source mappings remain separately gated. These
controls do not qualify I21–I24 or their predecessors or enable physics.

The compensated optical clock can retain more significant bits than a single
`long double`. The bridge now uses exact binary components and arbitrary-size
rational location, preserving the positive tail at the former projection
counterexample. The legacy `ExactPendingDepth`/scalar integrated API retain
their explicit projection refusal. Detailed repeated-flight, full-exponent
and independent Fraction controls are in
[`compensated_depth`](../compensated_depth/README.md). The
runner requires all 25 frozen named controls and reconciled counts/exits;
a zero-check or incomplete successful transcript is infrastructure ERROR.
