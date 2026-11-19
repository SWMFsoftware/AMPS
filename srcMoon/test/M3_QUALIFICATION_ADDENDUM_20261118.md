# Approved M3 statistical and numerical addendum

The user approved this design before new physical-cohort execution on
2026-11-18. [Machine-readable approved choices](M3_QUALIFICATION_ADDENDUM_20261118.json)
are SHA256-bound in the maintained scorer. The original
`M3_QUALIFICATION_PROTOCOL.json` remains byte-for-byte unchanged, SHA256
`29e183978e4c336891ea957726ddeb32e62898e7cb9eab8181909910b520ffaa`.
Old attempted results retain their original protocol and interpretation.

The numerical bias allocation is10% of the original ensemble gate: relative
0.002 for I21/I22 and0.003 for I23/I24. Original final2%/3%, direct source-rate
and stricter local kernel criteria are unchanged. Numerical verification must
use independently derived references at at least three predeclared steps and
one physical horizon; decreasing errors alone do not establish the allocation.
Quoted high-precision reference values are subtracted with Decimal before
normalization, avoiding a false zero caused by converting both operands first.

Original-parent endpoint survival is the primary homogeneous estimator.
I22 uses total-removal `lambda_hat=-log(S/N)/T`; I21/I23/I24 use its inverse
lifetime. Conditional source-channel counts give branch frequencies and partial
rate diagnostics. An event-time `D/sum(min(t_event,T))` estimate is retained as
a separately named diagnostic and never replaces a failing fixed primary.
Exact zero-hazard cases have a separate no-removal identity, with no invented
finite lifetime. Homogeneous cohorts refuse silent escape/surface censoring
or genuinely unequal initial trial weights. Children/splits are not independent
trials; physical packet/reference-weight tests remain separate numerical gates.

Linked branching uses simultaneous95% two-sided exact Clopper–Pearson intervals
with Bonferroni allocation over the ENTIRE frozen family. Every positive source
branch also needs interval half-width no more than10% of its independently
declared source probability. Source-zero branches require zero events and have
no invented relative error. U15's99% and U17's100000-event rule remain intact.
Source branches are not dropped, and surviving fractions are not renormalized.
The binary64 source-probability sum may differ from1 by its derived per-ratio
roundoff envelope; the planner records the original values/sum and performs NO
renormalization. Material differences are refused.

The fixed-cohort planner uses actual represented `rate*T`, not assumed exact
log2. Rate/lifetime gates have their own correct accepted-survival intervals.
To leave room for approved relative numerical bias `b` inside final gate `d`,
its stochastic planning margin is `(d-b)/(1+b)`; this does not tighten the
final acceptance. Exact binomial envelopes bound total removals and original
parent branch counts. Opposite monotonic Clopper–Pearson corners bound interval
width for every count in the envelope; complement symmetry protects efficient
planning for dominant branches. The family1% planning risk is allocated over
all endpoint, removal-envelope and branch-envelope conditions. Doubling yields
a sufficient fixed N, not a claimed minimum. All resource/search ceilings and
complete comparison-family sizes are explicit caller inputs. A cohort exceeding
them is deferred to an explicitly requested cluster package, never shortened.
The conditional branch-width planning confidence and the95% statistical
consistency test are distinct; failed seeds/intervals cannot be discarded.

`QualificationRandom` gives explicit parent/cohort/genealogy/purpose/counter
streams. Matched timestep/MPI comparisons share the intended physical keys;
rank/address/allocation/retry order never supplies seeds. The finite52-bit
midpoint grid, omitted end-support mass, optical-depth range and maximum draw
family are reported explicitly. The tail-support bound is NOT a universal
trajectory/model-error bound or mathematical independence proof.

See [qualification/CONTRACT.md](qualification/CONTRACT.md) and
[README.md](qualification/README.md) for exact entry points/commands/artifacts.
The readiness runner authenticates the plan/profile and any supplied actual
native receipt's eight artifact roles/token/layout. No current source-complete
physical cohort producer is admitted; a JSON Boolean, arbitrary CSV or all-green
administrative predecessors cannot create one. Formal I21–I24 still require
actual provider/source/numerical/ledger evidence and same-campaign predecessors.
This implementation creates no observation, native result or formal PASS.
