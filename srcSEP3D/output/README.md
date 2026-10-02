# srcSEP3D output interfaces


## Native shock telemetry

`shock_history.h` is a PIC/MPI-free stream boundary. `main_lib.cpp` captures
installed canonical SWCME state at a joined native clock boundary, compares
rank-min/max geometry/clock/speed and rank identities, and reduces owned
particle populations plus cumulative actual source-ledger allocations.
Reduced geometry is identical on every rank so the radius-stop branch is
collective even for a target near roundoff. Ghost particle copies are excluded.
`main.cpp` alone owns the rank-zero writer, broadcasts each I/O status before
another step, and closes before publishing `native-runtime.json`.

CSV starts at tick zero and contains every completed native tick. It records
physical `shock_active`; an inactive front is not fabricated into an active
shock. The writer rejects missing ticks, stale generations, changed provider
identity, clock mismatch, nonzero particles/injections and excessive rank
spreads. Existing histories are not overwritten. Interrupted streams may help
diagnosis but lack a completed runtime record and are not campaign evidence.

The runtime record contains geometry/clock/identity/completion facts only.
The Python validation runner owns input/executable hashes, launcher argv and
raw logs and adds a native evidence manifest only with independently frozen
absolute launch parameters. An unfitted control uses a relative-time PNG/EPS
and remains observationally unqualified. Propagation restart continuation is
currently rejected to avoid incomplete cumulative telemetry.

## SSE geometry and restart schema 4

New checkpoints serialize front kind, normalized propagation axis and angular
half width alongside apex kinematics. The little-endian payload schema is 4;
the file-family magic remains unchanged. Schema-3 checkpoints are still read
with Sphere geometry because older applications could only publish spheres.
Unknown geometry, invalid SSE axes/widths, checksum and identity mismatches
reject transactionally. Existing code/configuration identity policies remain
in force; accepting the old format does not bypass compatibility checks.

The existing propagation CSV columns `shock_radius_m` and `shock_speed_m_s`
mean apex distance/speed for SSE. They do not describe a radius common to all
directions or prove observer arrival. Use the frozen manifest's direction and
width and the directional surface geometry for flank arrival. Native MPI
telemetry checks geometry agreement as well as the epoch and apex state.
