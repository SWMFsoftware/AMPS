# Supporting normal environment/coefficient stage

ID: supporting extension of the existing M3 production-wiring owners; no new U/I ID.
Purpose: prove the extracted normal pre-particle stage retains frozen/evolving
SPICE epochs/SI vectors/frame coefficients and normal side effects; make the
same stage usable by finite cohorts without injection/sampling.
Entry: Moon::PrepareProductionParticleEnvironment, called by amps_time_step
between existing common Begin and PIC::TimeStep, and conditional native owner
between its common Begin and ExecutionTrackDefault. The complete original
orbit/thermal/cache/checkpoint block is copied byte-for-byte into the helper;
no second force implementation exists.
Configuration: normal Na/SPICE producer, frozen first/second physical-time
tracers at1/8MPI, plus legacy evolving mode at1/8MPI. Runtime inputs remain
explicit; actual global/block dt=.25s is an analytical clock-control fixture.
Evolving new tracers remain unsupported; use the maintained legacy mode to test
its unchanged evolving branch. Three stage calls/empty common steps per job.
Inputs: normal parsed application, same furnished kernels, actual initialized
mesh/material data; no observations. Actual M2 state stays disabled in these
controls; existing thermal side-effect source block is byte-preserved.
Oracle: independent furnished CSPICE SUN/EARTH/MOON state queries at declared
begin/middle/end epochs; Rdot*R^T angular velocity; actual body state matrices.
Units: SI m,m/s,rad/s; SPICE km and km/s converted once by1000. Time TDB ET.
Stochastic: none. No particles or random draws are needed by the stage fixture.
Procedure: normal amps_init, explicit consistent dt; frozen common Begin,
prepare, empty common Commit; legacy evolving prepare only. Repeat3times.
Metrics: complete actual vectors, angular velocity, frame blocks, source/start/
end epochs, sampling/particle counters; actual normal Na total acceleration
finite/nonzero on an allocated owner after preparation. Separate force/mover
fixtures retain integration/conservation scope.
Acceptance: exact same toolkit state values and declared frozen/end epochs;
rotation-derivative contraction agrees within256binary64eps*omega norm (a
predeclared arithmetic envelope for matrix contractions, never physical data
uncertainty); body derivative zero only for explicit frozen new-tracer policy;
no sampling/particle publication. No gate is changed after a failure.
Artifacts: source/build/executable identity; compiler and native logs; actual
furnished kernel paths/hashes; checks/metrics CSV; comparisonPNG/EPS; resultJSON.
Status: PASS means those software/wiring checks evaluated/met; FAIL evaluated
criterion missed; ERROR build/parser/runtime/provenance/resource failure;
SKIPPED missing normal producer/toolkit/capability. No formal M3 qualification.
Predecessor: fresh same-package normal M0 I01–10/nightly zero diff before runtime.
Resource: one original<=3600s package including code/setup/build/report/restore,
200GB process-tree cap,1/8MPI OMP1; original deadline is never reset.
