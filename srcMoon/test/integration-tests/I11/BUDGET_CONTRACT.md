# I11 production event budget and diagnostic-duration contract

This is supporting coverage for existing I11 and the existing U14/I02 Na
chemistry dispatch meaning. It assigns no new U/I ID, changes no physical
parameter or acceptance threshold, and cannot qualify I19 or M1 by itself.

| Field | Contract |
|---|---|
| Purpose | Explain actual neutral-Na gas accumulation independently of cell sampling and output writers, and identify the observables preventing stationary statistics. |
| Production entry points | Exosphere accepted source allocation; production accommodation callback; selected first/second-order AMPS mover escape and numerical-removal branches; actual configured Na reaction processor; replicated Exosphere prescribed-reservoir update. |
| Configuration | Approved 480 km frozen-epoch spherical M1 input, neutral Na, global species weight/time step, OpenMP=1, explicit MPI ranks and seed 731993. The deterministic linked fixture temporarily selects a lifetime-based chemistry step and controlled particle initial conditions, then restores the normal time step before transport. |
| Inputs / independent oracle | Exact model-particle entry-minus-exit identity and independent sum over allocated owned interior-cell particle lists. Weighted identity uses actual particle weights, including individual corrections; no theoretical source rate substitutes for an accepted event. |
| Units / frame / time | Gas-event counts are model particles; weighted inventory and events are neutral-Na atoms. Reservoir-supply count is the number of nonzero facet-prescription updates, not model particles; its weight is signed supplied atoms. Positions are LSO metres; frozen external epoch is TDB ET. Elapsed particle time is reported separately. Reservoir resets are informational external supply, never gas births. |
| Production scope | Emission by registered source ID, net surface adsorption, immediate reflection, photochemical parent loss/daughter creation, outer-domain escape, and known numerical mover removals. Missing paths produce an unexplained closure residual. MPI migration is neither birth nor physical loss. Full multispecies, evolving terrain and I19 closure remain separate. |
| Lifecycle | Enable a null-by-default, read-only production event observer. Event counters are cumulative and independent of visualization flags. Sampling occurs before injection/chemistry/motion in PIC::TimeStep: ledger rows record completed steps at that checkpoint. A final post-step snapshot covers events after the last sample. |
| Procedure | Complete one ordinary initialization step to populate the chemistry walker and frozen ephemeris; clear its particles before ledger accounting. Then create 4096 Na at a sunlit position and 4096 in cylindrical lunar shadow, allocated only by their MPI owner. Invoke the configured production chemistry dispatcher at dt=tau*ln(2). Check sunlit survival against p=1/2 and exact shadow survival. Then initialize transport fixtures (64 night-side impacts and one outward outer-boundary crossing), run 20 ordinary production steps with sources enabled, and check event closure. Run the diagnostic pilot only after these checks pass. |
| Chemistry oracle / criterion | Shadow survivors exactly 4096. Sunlit survival uses the two-sided Hoeffding 95% bound sqrt(log(2/0.05)/(2*4096)) around 1/2, declared before execution. This statistical fixture is not observational evidence. Neutral-only chemistry removes ionized Na; enabling a charged production case requires separately qualified product weights and is rejected by this Na-only adapter. |
| Budget criteria | Every checkpoint has exactly zero integer model-count residual. Weighted residual is bounded by gamma_n*sum(abs(terms)), gamma_n=n*u/(1-n*u), with long-double epsilon u and n conservatively including accumulated event additions, inventory additions and MPI reduction additions. Numerical removals remain separately visible and block treating the run as physically qualified. No empirical compensation is applied. |
| Observables | Per-metric raw history, lag-one correlation and its unchanged bound, three equilibrium changes, five-batch Student-t CI, and first failed gate. Postprocessing may reproduce the declared algorithm; it cannot discard a selected startup interval or change formal scores. |
| Pilot | Coarsest 480 km only, eight MPI ranks, 500-step samples, restart writes off; explicit wall budget (default 3600 s), 200 GiB RSS, 8 GiB host-available floor and 100 GiB free-disk floor. Wall-budget completion preserves a non-qualifying PARTIAL pilot with formal I11 NOT EVALUATED; resource/runtime failure remains ERROR. A candidate duration is never installed automatically. |
| Artifacts | Input and configuration snapshot, source/build/executable/SPICE hashes, stdout/stderr, ledger CSV, source-ID legend, per-metric CSV, statistical and budget PNG/EPS, result/manifest and process resources. |
| Expected PASS / failure modes | PASS is restricted to evaluated storage/chemistry/ledger fixture identities. Missing or corrupt evidence is ERROR; unavailable required capability is SKIPPED; an evaluated survival or closure violation is FAIL. A completed or partial pilot is never formal I11 PASS. |
| Data / predecessors | Configured hashed DE421 SPICE files; no observational package. Build, sampling-window prerequisite, linked budget/chemistry fixture, U07 and I10 precede interpretation of pilot durations. |

The November 12 source audit found that the configured Exosphere generic
reaction processor calls its own `TotalLifeTime` table, which has no Na
entry, rather than the configured Moon lifetime. The new linked fixture must
expose this discrepancy. Its production correction must retain the existing
Moon lifetime/shadow law and neutral-only species choice. Any resulting change
to the reviewed nightly regression is preserved for explicit reference review;
the checker must not replace its own frozen reference.
