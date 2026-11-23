# Stage 13 native material lifecycle: current boundary and reusable machinery

This audit describes the approved Sun-centred spherical `SHOCK_ONLY` first
profile and the actual checkout. It does not change admission. The existing
manufactured dipole/isothermal fixture verifies software and kinematics;
it does not supply a new wind/thermal closure or observational forcing.
The approved physical particle band is1.05 solar radii–5AU in HCI, with separate
numerical provider guards. Particle transport remains plasma-comoving and
uses one unchanged configured physical AMPS timestep between sampling outputs.

## Current implemented boundary

`adapters/swcme3d_spherical_material_adapter.cpp` captures an immutable actual
`CoronaSwcme::SphericalShockOnly::Event` as a neutral material interval. It
exports smooth upstream primitive/EOS/species/gradient data, diagnostic motion
and an independently evaluated local shock limit. It does not install PIC
vertex arrays, register a native line, own a source or publish a snapshot.
A downstream RH primitive supplies total pressure only: species heating,
downstream gradients and a downstream material volume are not inferred.
A sub-fast geometric surface has no fabricated jump or source eligibility.

The shared `Swcme3dSpherePfssParkerKinematicV1` and
`Swcme3dSpherePfssParkerThermalV1` profiles still require qualification for
`MaterialTransport` and `PhysicalDownstream` in
`src/models/sep_common/sep_material_background.h`. `RequireStep` and
`BindSnapshot` therefore reject actual live-step attachment. `S13AD01` exercises
both profiles using `InitializationLevel::None`. It provides linked adapter
evidence. The thermal profile supplies an independent identity and the actual
electron, proton and alpha temperatures and Cartesian gradients at query time.

The production `amps_time_step()` path in `main_lib.cpp` currently performs its
material geometry transaction only through `ReducedShock::Native::Advance`.
The spherical Stage13 validation selector reads named cases and checks input
readiness; it does not secretly construct a different native provider.
Neither a successful core trajectory nor empty particle storage is a reason
to bypass the missing physical closure or this native ownership boundary.

## Reusable zero-particle geometry path after a wind realization is qualified

The reusable container is `PIC::FieldLine::cFieldLine`, not a new Cartesian
mesh. The following mechanisms already exist and should be reused:

1. `Geometry::Prepare` in `src/pic/pic_field_line_geometry.h` accepts a complete
   position-map callback and an inner magnetic tracing callback. It moves
   every supplied old vertex, prepends material vertices and retires only
   wholly external segments using the actual chord minimum radius.
   `Geometry::TraceAlongField` traces the target-epoch B direction; material
   ordering resolves the sign independently of magnetic polarity. This does
   not equate solar-wind motion with motion parallel to B.
2. `cFieldLine::RegisterGeometryOwner` records persistent physical identity,
   epoch and material labels before sources/sampling start. Preparation checks
   that identity, epoch and generation. It does not fabricate their values.
3. `cFieldLine::PrepareGeometryStep` in `src/pic/pic_field_line.cpp` stages
   vertex and segment data, index/ownership resources, sample retirement and
   particle remapping before changing the live line. New vertex history must
   come from the actual provider's backward characteristic and departure
   state, or be explicitly unavailable. Copying a neighbour's past is not
   initialization of a new material record. New simulation segments start
   empty and use the existing source rules.
4. `CommitGeometryStep` rejects a foreign/stale/consumed plan, altered departure
   coordinates, changed particle lists/weights/state, changed sampling exposure
   or stale ownership resources. The segment `SamplingExposure` counters stay
   attached to their actual material segment. Empty segments still accumulate
   one sampling opportunity per full physical step.
5. `SnapshotStore::PreparePublication/CommitPublication` preserves the existing
   single publication authority and read-phase exclusion. The collective
   preparation/coordinate-agreement pattern in `reduced_shock_geometry.cpp`
   demonstrates how to stage candidate topology, primitive data, connectivity
   and metadata before a joined MPI commit. This is a pattern to reuse, not an
   authority to label a spherical owner `ReducedShock` to satisfy its guard.

A future background-only Stage13 driver must first hold the actual qualified
immutable wind/event, build its real B-traced line, and stage all current and
previous primitives at the corresponding physical epochs. It must use the
real provider topology/sector at every trajectory and derivative query,
including guard and inner-entry paths. It must preflight the complete interval,
then commit geometry and matching snapshot/connection generations after all
workers and readers have joined. No particle source, wave update, RNG event,
clock mutation or provider switch is implied by a zero-particle geometry test.
The global allocated-particle count must remain exactly zero on every accepted
background-only step, including newly inserted segments.

There is an explicit fixed-duration dependency. The core offers
`Event::AdvanceFixed(x,begin,end,dt,controls)` to integrate the supplied physical
duration, while its diagnostic `Move(begin,end)` uses `end-begin`.
The detached neutral owner currently uses diagnostic `Move`; at a large
absolute epoch that difference can differ from the intended dt. The future
native geometry callback must capture the authoritative frozen driver-window
duration and use `AdvanceFixed` for the complete host step. `Geometry::Context`
currently carries endpoint epochs rather than a physical-dt authority.
No shortened step, endpoint-difference substitute or sampling-dt change is
permitted. Geometry integration subdivisions are numerical work within the
same fixed physical step, not extra physical sampling opportunities.

## Missing native qualification before live transport/source admission

- **Physical wind and thermal realization:** qualify the actual event's
  nonradial open tubes. Controlled tests cover continuity, signed induction,
  the field-aligned momentum equation and species heat accounting over a
  resolved axisymmetric band. The full transverse force and observational
  thermal/work bounds still require event inputs.
- **Native startup/publication:** install one real spherical material owner,
  its actual tube metrics, all mandatory vertex datums and immutable physical
  identity/configuration. Agree preparation across ranks before publication.
  A descriptor or a renamed reduced event cannot substitute for the owner.
- **Source/connectivity:** reuse stable source-line namespaces and keyed
  counters, but query the actual spherical front and moving line at the real
  birth epoch. Preserve the already reviewed spectrum/frame/normalization and
  shock-return policy. Validate accepted fast-shock classification, one-sided
  state, disconnection, source area/measure and complete signed SI ledgers.
  Local RH pressure does not select a downstream species-temperature law.
- **Physical escape:** implement/check first passage through the physical5AU
  sphere while the provider retains its larger numerical guard. Retiring an
  entire outside segment does not by itself prove first-passage escape for a
  particle inside a partially intersected segment. The larger-domain control
  must retain common source support and the same physical dt.
- **Particle remap and sampling:** check field-line coordinate ownership,
  derived Cartesian projection, tube-volume measure, inserted empty material,
  retired segment numerators/exposure and all-rank conservation using the
  actual new provider. Preserve the three existing particle mover equations.
- **Restart:** the maintained material checkpoint reconstruction currently
  owns a reduced-core interval and source ledger; fixed-analytic checkpoints
  do not admit an arbitrary moving material owner. Add the actual spherical
  event/wind identity, geometry/material labels/generations, primitive history,
  sampling exposure, source counters and boundary/ledger state to a reviewed
  continuation path. Reject unsupported/foreign/stale state before mutation.
- **Observational readiness:** retain authentic instrument/source uncertainties
  and explicit missing inputs in the Stage13 case manifests. An OV01
  continuity PASS, synthetic provider PASS or source checksum is not a
  qualification of new spherical SEP transport or all SWCME3D shapes.

The roadmap's complete Stage13 exit additionally includes nonspherical normals,
speeds, areas, finite support, tangencies and material budgets for each future
shape/profile. The approved spherical first profile does not admit ellipsoid,
finite-SSE, downstream volume or FULL_ICME by inference.

## Focused implementation in this audit

The existing Routine `S13AD01` callback was expanded in
`test/stage13_material_adapter_validation.cpp` without changing its ID,
initialization level, input fixture or admission policy. It checks exact bridge
translation of U/B, pressure/EOS, all configured species and both tensors;
independent simultaneous capture generations; opposite declared sectors;
retained immutable-owner lifetime; actual snapshot-binding/downstream rejection;
wrong branch/side/time and absent capability rejection; diagnostic motion
parity; and genuine sub-fast classification without an invented downstream
state. Exact equality applies to a value bridge, not an observational tolerance.
Its assertion count and evidence scope are reported in normal native results.

## Regression commands and scope

After rebuilding a matching native executable, run the actual new callback:

```sh
mpiexec -n 1 /absolute/current/amps --test S13AD01 \
  --test-json /absolute/fresh/S13AD01.json \
  --test-junit /absolute/fresh/S13AD01.xml
mpiexec -n 2 /absolute/current/amps --test S13AD01 \
  --test-json /absolute/fresh/S13AD01-mpi2.json \
  --test-junit /absolute/fresh/S13AD01-mpi2.xml
make -C src/models/sep_corona_swcme -j4 test-spherical-shock-only
make -C srcSEP test-stage13-cases SEP_EXECUTABLE=/absolute/current/amps
make -C srcSEP test-material-background-unit
make -C srcSEP test-material-tube-kinematics-unit
```

Those commands verify adapter/core/input and retained shared software contracts;
they do not qualify a missing native Stage13 lifecycle. When implementing that
lifecycle, retain the existing native geometry/source/mover and checkpoint
protocols as compatibility gates:

```sh
python3 -B srcSEP/test/run_native_standalone_integration.py \
  --build-root /absolute/current/native-build --ranks 1 2 4 8 \
  --output-dir /absolute/fresh/native-prepared --prepare-only
python3 -B srcSEP/test/run_material_checkpoint_tests.py \
  --native-prepared /absolute/fresh/native-prepared \
  --output /absolute/fresh/material-checkpoint --prepare-only
make -C srcSEP -j1 test test-controlled-analytical \
  SEP_EXECUTABLE=/absolute/current/amps
```

`--prepare-only` builds and emits launch manifests; it is not a test PASS.
Execute every declared native job with its original seed, rank/thread settings
and fixed dt, retain actual process receipts, and run the corresponding
`--collect` stage only after completion. New Stage13-specific geometry/source,
escape and continuation tests must consume the new owner; passing the retained
reduced-owner tests alone cannot admit the new profile.
