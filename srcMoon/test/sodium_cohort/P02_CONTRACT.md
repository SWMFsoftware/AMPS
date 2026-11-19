# Independent P02 mechanical-error verifier — contract before implementation

Local supporting verification under the existing U15/I21 sodium scope; no new
permanent ID and no observational or whole-M3 qualification claim. The approved
budget is an explicit path/SHA descriptor to `M3_P02_NUMERICAL_BUDGET.json`.
Its four values are read, never duplicated as numerical defaults in code.

Inputs are a frozen JSON binding of an authenticated source profile, approved
budget and original completed native wrapper owner manifest/result descriptors.
Each descriptor carries an exact path/SHA. Native input receipt, complete input
include graph, parent CSV and every rank's hexadecimal event file must match the
original wrapper artifact map. Missing data/capability is SKIPPED, corruption or
schema/hash/owner contradictions ERROR, evaluated bounds exceeded FAIL. An
already evaluated failure cannot disappear behind a later missing input.

Production interpretation is `ReducedKinematics::EvaluateModelErrors`: source
reference mass m_r is the sum of neutral isotope masses minus integer charge
in electron masses, physical mass m_p is the independently sourced state CSV
value. Relative inertia error is |m_r−m_p|/m_p, relative acceleration error is
|m_r−m_p|/m_r. These denominators must not be interchanged. Kinetic discrepancy
uses the SAME mass in K_class=m_r*v²/2 and K_rel=(gamma−1)m_r*c²; separate mass
error is not folded into its denominator. Thus its relative value is
|K_rel−K_class|/K_rel, with the exact zero-speed limit0. SI c is authenticated.

The implementation will independently evaluate the kinetic inequality using
rational beta² and a squared nonnegative inequality, avoiding underflow or
subtractive sqrt cancellation. Decimal arithmetic is diagnostic only. No
production math, model-error routine or verified boolean is imported.

Actual parent/electron inertial velocities are taken from exact native hex
records. The actual committed ion velocity is recorded in LSO; its inertial
velocity is reconstructed by EXACT rational inversion of the recorded affine
R and offset, not an assumed zero velocity or a guessed frame. Matrix
orthogonality/proper determinant must satisfy the retained source rotation
bounds. The authenticated original input must name SUN-origin inertial axes;
this verifier does not certify SPICE frame identity or create a new physical
frame gate. It checks the explicit speed cap for both recorded event velocities
and admitted inertial velocities, while the relativistic discrepancy is labelled
for the declared inertial frame only. It checks the worst-case declared speed
cap independently as well, matching production's conservative admission.

Scope is specifically the approved 23Na neutral ground, 23Na+ ground, free
electron channel and source-effective-photon policy. Source profile, state IDs,
physical mass source roles, isotope constants and complete input graph remain
hash-bound. The declared native state-specific caps must equal the pinned
approved budget values; no arbitrary receipt bound supplies its own approval.

Tests: labelled pure rational/Decimal arithmetic fixtures for differing mass
error denominators, zero and near-boundary kinetic cases, subnormal speeds,
limits and invalid numbers; exact affine inverse and proper-rotation controls;
missing/corrupt/changed source or original-owner administrative binding controls.
Real retained profile controls have `requires_authoritative_source=True` and
remain source QA, never fabricated native event/observation evidence. An actual
physical owner is needed for an actual event-scope result. No fixture is written
as a real producer receipt or treated as a physical qualification.

Outputs in a NEW directory: per-criterion result.json, comparison.csv,
run_manifest.json and scalar check log. Exact fractions preserve binary hex and
source decimals, with explicit source/budget/verifier hashes and native owner
paths. No plot is needed for scalar guard checks. Empty actual event sets retain
source checks but cannot be called observed event-speed verification. No files
in canonical validation data or the shared installed database are modified.

The first actual43-event audit exposed a verifier implementation error: it
scaled off-diagonal Gram residuals by the near-zero dot product instead of
the dimensionless unit-vector norm. Production `ValidateFrame` declares the
column Gram test against identity with scale1, and its determinant test also
uses scale1. The correction uses that SAME unit scale and unchanged native
absolute/relative bounds. Tests freeze the actual rounded matrix and exact
1e−15 off-diagonal boundary plus a strictly larger rational violation. The
initial failed audit is retained; no event, source or threshold is changed.
