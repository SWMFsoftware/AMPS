# Compensated variable-history optical depth

Written before the production correction, 2026-11-18. This extends the existing
U15/U16 numerical mechanics; no permanent U/I ID or physical threshold changes.

The exact two-binary64 pending depth can require over1000 exponent positions
(for example1 minus the smallest positive double). The existing preparatory
bridge refused projection to long double. Wiring that preparatory limit into
the real variable controller caused ordinary repeated7e-6s^-1 flights to stop
after their first no-event prefix. Preserve the old explicit projection API,
but use a genuine full-component interface for actual variable histories.

Production entry points: `PiecewiseHazard::LocateEvent` (new compensated-depth
overload), `ReactionExposure::{Locate,Evaluate}`, the SAME instantaneous
`ReactionExecution::SelectHazardChannel`, and an explicit located-result clock
publication adapter. Existing scalar and constant-history interfaces retain
their meanings and numerical paths. SamePendingDepth compares both serialized
pending components for the parent's disposable numerical cache.

Every finite binary input is converted exactly to an arbitrary-size rational
using the existing production arithmetic. Comparison and first-hitting-time
certification never use a rounded scalar sum, fixed256-bit mantissa, a CDF
tolerance or a sign/zero repair. The entire history support is checked before
any event decision. Exact segment-end events retain left-limit semantics.

No-event residuals are calculated from the exact rational pending depth minus
the exact rational affine integral. Native storage remains the existing two
binary64 components. Round the high component, then its exact signed residual;
retain a rigorous explicit bound for any remaining storage rounding. This is
a numerical representation receipt, NOT a physical acceptance allowance or
permission to move an event. Zero rounding bound means exact representation.
Positive residual underflow/overflow or invalid state is refused transactionally
instead of turning it into an event or dropping a low component. Event flags
come from the exact locator, never comparison with a displayed rounded area.

## Test specification

| Field | Specification |
|---|---|
| Purpose | Actual full-component source-clock continuation and exact near-boundary first crossing |
| Configuration | Pure C++17 production kernels, no AMPS/PIC/MPI build or data files |
| Inputs | Predeclared mathematical clocks, affine rates and intervals; no observations |
| Independent oracle | Python Fraction on exact binary components and closed-form/rational affine integrals |
| Units/frame/time | Dimensionless depth, rates s^-1, relative physical seconds; no lunar frame/epoch selected |
| Stochastic controls | Explicit open uniforms; no RNG or redrawing |
| Procedure | Reproduce original failure;10000rate7e-6 flights; equivalent dyadic partitions/outer intervals; positive/negative low components; full double exponent-span case; non-dyadic residual storage bound; invalid/underflow refusal |
| Metrics/criteria | Exact rational event/no-event decisions; first representable non-early time; exact conservation of dyadic cumulative depths where two-component representable; non-dyadic residual error bounded by recorded storage rounding; unchanged outputs on refusal |
| Expected PASS | Every evaluated arithmetic/provenance/transactional criterion met |
| Expected failures | Discarded low component, rounded early event, hidden truncation, wrong instantaneous CDF, state mutation on refusal |
| Required data/hashes | No science data; every production/probe/script/compiler/executable identity retained |
| Predecessors | Source/local mechanics only; formal M1/M2 and physical I21–I24 gates unchanged |
| Artifacts | Native stdout/stderr, compiler log, result/run manifest, Fraction comparison CSV and PNG/EPS diagnostics |
| Status | PASS evaluated/met; FAIL numerical criterion missed; ERROR build/parser/runtime/infrastructure; SKIPPED absent compiler/plot dependency |

Qualification must still account for accumulated numerical error by the
unchanged convergence protocol. This correction supplies neither a physical
reaction provider nor source/observational readiness.
