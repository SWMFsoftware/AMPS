# Supporting native Moon force/tracer integration contract

Pre-implementation contract, 2026-11-17. Existing U03–U06 force oracles and
M3 ion-mover/integration coverage are extended without a new permanent ID.
**This is frozen-snapshot production wiring/integration verification, not
independent physical force validation, observational validation, a measured
plasma driver, or completion of formal I21–I24.**

| Field | Frozen supporting contract |
|---|---|
| Purpose/entry | Both normal `Moon::MoveParticle` order selectors must execute actual `Moon::TotalParticleAcceleration` through core ForceCutSurface's NULL acceleration callback. No custom force is installed into the mover. |
| Production build | Same successful clean eight-species Moon library and generated ABI as native particle-backend fixture; individual weights/local species timesteps/local weights ON, unrelated velocity sampler OFF, MPI-only/OMP1. Retain producer library/build receipt and test-only executable identity separately. |
| Inputs | Analytical, observation-free frozen LSO snapshot: Moon→Sun1AU along+X; Moon→Earth(−3.8e8,2e7,5e6)m; frozen object velocity(0,30000,0)m/s; frame angular velocity(1e−7,−3e−7,2e−6)rad/s. Prescribed E=(1e−6,−2e−6,3e−6)V/m, B=(3e−9,−4e−9,5e−9)T, explicit frameLSO. These are dimensioned numerical controls, not claimed ephemerides or observations. |
| Cohort | Exactly one neutralNa and oneNa+ global tracer, with the same initial x=(2.11e6,123456,234567)m, v=(730,−410,95)m/s and represented number1. No reactions/collisions/surface impacts/injection. A64-leaf real distributed mesh contains the complete path; each particle starts only on its actual MPI owner. No MPI migration claim follows when the path stays in its initial leaf. |
| Frame/time | Frozen stationary LSO physical fields and force globals; arbitrary declared initial physical epoch1000s, requested local interval4s. No UTC/TDB/actual-epoch interpretation. Node species timesteps intentionally differ from4s; caller interval controls advancement. |
| References | At initial state independently form lunar gravity, differential Sun/Earth gravity, centrifugal/Coriolis terms and q/m(E+v×B); radiation wiring uses the existing compiled Combi scalar table and independent direction/velocity projection. Independent long-double RK4 integrates the **actual production total acceleration**; it is an integrator reference, not an independent force law. U03–U06 remain required for independent term qualification. |
| Numerical design | For each species/order use the same T=4s and fixed maximum substeps h=1,.5,.25s, adaptive OFF. Reference RK4 uses4096 and8192 uniform steps, long-double state accumulation and actual double production-force calls. Selection frozen before outputs. |
| Metrics/acceptance | Every declared force contribution is finite/nonzero at input; total initial production force matches the independently assembled vector within5e−12 relative or1e−13m/s² absolute (rounding-level wiring bound). Separate position[m]/velocity[m/s] error norms decrease; observed orders from consecutive refinements must be[.8,1.2] for first order and[1.8,2.2] for second. Reference successive difference must be≤0.001 times finest tracer error, separately for position/velocity; otherwise ERROR/insufficient oracle resolution, not relaxed acceptance. |
| Invariants | Actual compiled species and represented N preserved; stored opaque incoming optical clock bytes unchanged with event callbacks absent; actual trace context start1000/request4/elapsed4/remaining0/order correct; stored physical epoch1004; global cohort count1 for each species. |
| Procedure/artifacts | Thin probe only against successful retained producer library; no application recompile. Execute order1/order2 at1/8MPI, preserve source/hash/compiler/library/generated input identities, stdout/stderr, source-defined control config, per-case state/errors/order CSVs, named checks/result/run manifests and diagnostic PNG/EPS comparisons. |
| Data/predecessors | No external observational package/kernel consumption. Actual successful backend producer build needed. Formal predecessor gates and source/data readiness remain separate. |
| Status | PASS declared supporting checks only; evaluated integrator/wiring miss FAIL; compile/launcher/partial evidence/reference-insufficiency failure ERROR; absent qualified retained build SKIPPED. Never promote shared-force RK4 agreement to physical validation. |

Snapshot vectors are explicit controls; they do not fit any residual. No
constant B table is modified. Both explicit tracer orders obtain immutable
fields through the ordinary Moon runtime field provider. The script retains
its test source separately because adding a thin probe does not change the
already compiled production object/library bytes.


The source boundary extension also calls both actual configured injection
macros on a freshly allocated particle with no epoch. A stationary local
species step of0.125s supplies remaining0.03125s; the independent physical
identity is birth plus remaining equals the declared step endpoint. Both
orders/layouts must retain the inferred birth1000.09375s, endpoint1000.125s
and selected order without a second offset draw. The binary fractions are
exact numerical controls, not measured lunar times. Existing force/order
criteria remain unchanged; the eight additional check evaluations are
reported separately from older81-check source snapshots.
