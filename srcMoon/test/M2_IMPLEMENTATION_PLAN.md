# M2 development before M1 closure

The user approved this development-order exception on 2026-11-12. It permits
implementation and supporting verification, not promotion of M2 or bypass of
I13. The source roadmap amendment controls this exception; the common roadmap
still controls physics, acceptance, data interpretation and formal milestones.

## Inventory and implementation contract

| Field | Contract |
|---|---|
| IDs / purpose | Supporting existing U08-U13 and I14-I20; no new permanent IDs. Complete missing M2 production accounting and strengthen analytical geometry/thermal checks before the expensive formal campaign. |
| Production capabilities | Actual Exosphere source allocation, surface adsorption/desorption and reservoir updates; AMPS mover losses; production LOLA/BVH/thermal services. |
| Entry points | `Exosphere::SourceProcesses`, Exosphere surface population/flux balance, `PIC::ParticleEvents`, `Moon::SurfaceTemperature`, `Moon::SurfaceEnvironment`, configured movers and ordinary `amps_time_step`. |
| Configuration | Dedicated neutral Na/Ar M2 build, photolytic chemistry OFF as already declared by I19; normal production callbacks, seed 731993, eight MPI ranks/OpenMP=1 for bounded supporting runs. Canonical neutral-Na baseline is restored after builds. |
| Inputs | Reviewed M2 application receipt and existing species parameter contracts; native D01/D02/D03 files only for explicitly data-dependent checks. No assumed observations or new physical coefficients. |
| Independent oracle | Signed gas and surface atom conservation, mass conversion using production species masses, exact source/sink transfer identities; analytical geometry limits and conduction invariants; native archive labels and frozen hashes for data QA. |
| Units | Gas event counts are model particles; weighted events/inventories are atoms; mass is kg; rates use synchronized simulation seconds. Coordinates are LSO or explicitly identified MOON_ME_DE421 metres; thermal quantities are K, W/m2 and J/m2. |
| Time / frame | SPICE TDB internally, UTC input; freeze/evolve exactly as stated by the archived receipt. No fabricated observation geometry. |
| Stochastic controls | Explicit seed, particle weights, ranks/threads and sample counts retained. Local analytical fixtures are verification inputs, never observations. |
| Procedure | Inventory paths first; instrument actual accepted events and actual reservoir changes without changing RNG or physical laws; compare independently summed owned gas and replicated logical-surface states; run existing local geometry/thermal gates; audit D02/D03 requirements; preserve separate formal predecessor-blocked results. |
| Metrics | Gas/surface/combined signed residuals by species, event counts/weights, source classification, numerical loss, analytic errors, data hashes/readiness and exact missing prerequisites. |
| Acceptance | Existing I19 deterministic precision and <0.5% dominant-source imbalance remain binding; no numerical loss hidden as escape. Local U08-U13 keep existing criteria. Supporting PASS applies only to the evaluated identity, never to formal I19/M2. |
| Expected results | Supporting checks PASS when their exact criteria are evaluated and met; formal I14-I20 remain SKIPPED without I13. Data packages retain real PARTIAL/INCOMPLETE states. |
| Failure modes | Missing event path; inconsistent weighting; adsorption counted twice; reservoir copied across MPI ranks; mismatched accounting times; physical/numerical loss conflated; malformed native labels or hashes; runtime/resource failure. |
| Required data / hashes | None for analytical bookkeeping; configured SPICE hashes for linked runs. Data checks read dataset README/provenance/QA and verify native raw/processed hashes before use. |
| Predecessors | M0 baseline and applicable local gates support development. Formal qualification additionally requires unchanged I11-I13 and dataset/operator gates. |
| Artifacts | Separate logs, input/source/build/executable/data identities, seeds/layout, row-level CSV, diagnostic PNG/EPS, supporting result/manifest, formal SKIPPED result, resource history and aggregate counts. Earlier results are never overwritten. |
| Status semantics | Supporting PASS/FAIL/ERROR/SKIPPED describe only declared supporting checks. A preflight PASS is not runtime PASS. Formal I-series status stays SKIPPED when its predecessor is unavailable. No observations/hold-out residuals are scored by this development campaign. |

## Work order

1. Extend the existing read-only event observer and connect actual Na/Ar gas
   and surface transfer paths. Preserve prescribed-reservoir external supply
   separately from finite flux-balanced inventory. Do not infer source events
   from theoretical rates or silently compensate for unexplained residuals.
2. Add a bounded supporting linked fixture and independent signed-budget
   evaluator. Exercise adsorption, desorption, external injection and escape;
   report any path that is unavailable instead of claiming coverage.
3. Run/extend existing U08-U13 analytical tests through production services.
   Do not duplicate the model's force, intersection or thermal laws in tests.
4. Add reproducible D02/D03 package audit/qualification support. Native
   metadata, collection-specific processing QA and frozen selections control
   readiness; wrappers alone cannot make a package READY. Preserve raw data
   and genuine gaps; never expose December residuals in development.
5. Register one bounded development campaign in the global runner, retaining
   all individual logs/results and supporting-vs-formal counts. Keep the
   existing qualification path and I13 prerequisite unchanged.

M1 closure still requires a reviewed spin-up/common-duration/statistical
protocol and I11-I13 PASS. I20 remains blocked by I15/I17/I18/I19; this work
does not claim to finish every M2 observational or numerical gate.
