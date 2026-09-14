# Step 14 change manifest: canonical source cleanup

Step 14 removes the legacy mover implementation surface and turns its absence
into a maintained acceptance condition.

## Removed

- monolithic `mover.cpp`, legacy `fte_mover.cpp`, and stale
  `fte_mover_dmumu.cpp`;
- deprecated mover declarations and the mutable `ParticleMoverPtr` dispatch
  escape hatch;
- all mover-name compatibility aliases;
- inactive transport-coefficient/event-scattering globals;
- generated binary, coverage, report, output, and nested-archive paths from
  future source packages through `.gitignore` and `DOC03` checks. Step 15
  extends this maintained boundary to Python bytecode caches and campaign
  report names.

## Resulting layout

- `parker_mover.cpp`: canonical pitch-angle-averaged Parker adapter;
- `focused_transport_dmumu.cpp`: canonical coefficient-driven FTE adapter;
- `focused_transport_mfp.cpp`: canonical event-driven MFP adapter;
- `mover_state.cpp`: only configuration shared by those three adapters;
- `production_mover_runtime.cpp`: private three-entry implementation mapping
  behind the validated dispatcher.

`MIGRATION_MANIFEST.md` maps every removed name and C++ symbol to an explicit
replacement or to the separate Cartesian application. `test/run_step14_tests.sh`
implements `DOC01`–`DOC03`, `WARN01`, and the source-only portion of `SAN01`.
