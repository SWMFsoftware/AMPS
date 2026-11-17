# Approved M2 execution and numerical protocol implementation

The user approved the v3 proposal and requested implementation on 2026-11-13.
This contract precedes production/test changes. The approved scientific error
budgets are numerical verification criteria, not observed uncertainties.
The external data tree remains read-only; no commit or push is authorized.

| Required field | Implementation contract |
|---|---|
| ID/name | Supporting I19/I20 protocol and existing M2 development phase; preserve U01-U26/I01-I46. |
| Purpose | Execute bounded production verification within one aggregate hour; implement independent controls/diagnostics and separately budgeted cluster qualification. |
| Production capability | Normal Na/Ar injection, moving LOLA boundary, sticking/release, force mover, finite surface stock, thermal state, particle/surface accounting. |
| Entry points | MoonIntegrationTests, MoonM2Budget, main_lib::amps_time_step, Exosphere finite injector/surface interaction, MoonInput parser, existing global/development/build runners. |
| Configuration | Dedicated Na/Ar, chemistry OFF, individual weights ON, eight MPI ranks/OpenMP=1; optional diagnostics disabled in ordinary M0/M1. |
| Inputs | Explicit .in numerical case/seed/dt/weight/cohort controls and approved protocol; immutable D01 and SPICE; D02/D03 preflight without residual inspection. |
| Independent oracle | Exact signed gas/surface identities, independent first-order hazard integrals, constant/piecewise warming limits, separately implemented confidence calculations and process/deadline protocol fixtures. |
| Units | SI seconds/TDB epochs, atom packet weights, physical atoms/kg, hazard s^-1, surface temperatures K; wall seconds versus allocated CPU-core hours explicit. |
| Frame/time | MOON_ME_DE421 surface, LSO particles, SPICE TDB; common Na/Ar interval and integer checkpoint indices; no invented D02 epoch. |
| Stochastic controls | Actual configured seeds/rank seeds and RNG recorded; no RNG in diagnostics; tag fields travel in real MPI particle records. |
| Procedure | Verify optional parser/observer/control defaults, implement all-process deadline, build privately, run local/linked cases and short two-seed 3x3 sweeps; preserve formal predecessor/data SKIPPED. Cluster jobs remain explicit and cost/gate checked. |
| Metrics | Real wall/build/plot/cleanup durations; per-case seed/dt/weight/epoch, exact counts and roundoff-bounded signed atoms, first-release/censored residence/occupancy, CPU-hours/resources, scoped statuses/artifacts. |
| Acceptance | 3600-s whole development envelope, configured cleanup reserve; exact/frozen conservation gates; approved 1% numerical/one-percentage-point fraction/1%-H residence and 0.1% long-case loss margins used only for the full declared qualification design. |
| Expected PASS | Evaluated supporting kernel/wiring/accounting checks pass; short two-seed sweep never substitutes for full stochastic or month-long qualification. |
| Failure modes | Orphan MPI/build process, reset child deadline, stale/reused seed, dt-dependent weight, duplicate replicated stock, dropped tag/mass, unsupported thermal/data state, corrupt inputs, numerical failure or exhausted resource. |
| Required data/hashes | Existing qualified D01 LDEM_4 and SPICE hash lists verified by existing manifests. Formal I20 additionally needs exact READY D02/D03 and metadata; missing D02 stays SKIPPED. |
| Predecessors | Local analytical gates and baseline regressions for development; same-campaign I13 and unchanged I15/I17/I18/I19/data gates for formal I20. |
| Artifacts | New campaign directory, stdout/stderr, result/summary/manifests, input/config/source/build/executable/data hashes, ledger/cohort CSV, CSV/PNG/EPS comparisons, seed/layout, resource and CPU-hour accounting, all attempts. |
| Status semantics | Selected evaluated miss FAIL; runtime/build/parser/deadline/resource fault ERROR; missing/unselected long qualification SKIPPED; supporting PASS cannot open formal gates. |

Pre-change source revision: 06eafde4802f0344694612c8f3c21f2b28e1a09f.
Source roadmap hash: bde5f001bdcbe072d7a8452d3f2afec93856f50cd1146560da3cd3027a96d2dc.
External roadmap hash: d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4.
Source amendments control approved development/execution; common science and
native metadata control qualification. Record subsequent amendments/hashes;
never rewrite the hashes in historical run evidence.
