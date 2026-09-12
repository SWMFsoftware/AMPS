# SWCME semi-analytical solar-wind / CME-shock model

SWCME provides lightweight 1-D and 3-D solar-wind/CME-shock backgrounds for
energetic-particle transport studies.  The 3-D model combines a Parker magnetic
field, Leblanc density profile, analytical CME/shock kinematics, configurable
shock geometry, and phenomenological sheath/ejecta fields.

## Canonical defaults and declared model scope

`swcme_defaults.hpp` is the single source of truth for dimensionality-independent
SWCME defaults.  A default-constructed 1-D and 3-D model now use the same
ambient plasma, Parker normalization, apex kinematics, region parameters,
smoothing widths, and source weighting.  In particular, the common baseline is

```text
V_sw                    = 400 km/s
n(1 AU)                 = 5 cm^-3
|B|(1 AU, reference)    = 5 nT
T_p                     = 1.2e5 K
gamma                    = 5/3
Parker reference sinθ    = 1 (equatorial normalization only)
Omega_sun                = 2.86533e-6 rad/s model convention
kinematics               = DBM
DBM reference radius     = 20 R_s
V0                       = 1500 km/s
Gamma                    = 1e-7 km^-1
regions                  = SHOCK_ONLY
shock acceleration       = SOURCE
```

The default 3-D science geometry is the finite `SSE` cap with 40-degree half
width. `Sphere` and `Ellipsoid` remain available, but a Sun-centered sphere is
primarily a verification/reduction geometry and is no longer the default science
front.

The Parker magnetic convention is explicit. `B1AU_nT` is a **positive total
field magnitude** at 1 AU at the documented reference latitude
`sin(theta_ref)=1`; the common solar-wind preparation converts that value to the
radial `Br(1 AU)` normalization.  Positive radial field means outward polarity.
In 3-D, the public `sin_theta` member remains only for source compatibility as
this normalization latitude; it does **not** control local Parker winding.  Local
`sin(theta)=|Omega_hat x e_r|` is always obtained from geometry.

Two model scopes are named explicitly:

- `CONTROLLED_SEP_PRE_SHOCK` = `SHOCK_ONLY + SOURCE`.  This is the canonical
  baseline for the SEP connectivity/perpendicular-diffusion study.  At a given
  observer, its Parker/Leblanc transport background is declared in scope only
  while the modeled shock has not reached that observer.
- `FULL_ICME_DIAGNOSTIC` = `FULL_ICME + RESOLVED_COMPRESSION`.  This is the
  optional phenomenological sheath/ejecta diagnostic model.  It is not a claim
  of a validated global ICME magnetic structure; for example, the ejecta field
  remains Parker-like.

The scope is *derived* from the validated region/acceleration pair rather than
stored as a third independent option.  `Model::model_scope()` therefore cannot
disagree with the actual transport/source representation.  Both dimensional
models expose `observer_scope_status(...)`; the 3-D implementation uses the same
finite shock geometry as the production shock/connectivity routines, so an
observer outside an SSE cap is not falsely marked post-shock merely because the
apex has passed its heliocentric radius.

`resolved_configuration_manifest(const Params&)` in both public namespaces emits
a deterministic complete key/value snapshot of the resolved configuration,
including inactive compatibility parameters, data-driven tables, geometry and
the derived model scope.  This function is intended as the configuration block
for the campaign/run manifests emitted by `test/run_tests.py`; event-specific
overrides should never exist only in ad-hoc driver code.  The schema/convention
version is `SWCME_CONFIG_VERSION = 2`.

## SWCME-to-SEP / AMPS source interface

`swcme_sep_source.hpp` and `swcme_sep_interface.hpp` are the production-facing
integration layer between the validated SWCME background/shock model and an SEP
transport code such as AMPS.  The adapters do not duplicate CME, Parker,
connectivity, or Rankine-Hugoniot physics: they consume the same prepared
`StepState`, checked background evaluators, local shock state, corrected shock
mesh, and cobpoint solver exercised by the validation suite.

### Prepared-step usage

Both dimensional adapters expose the same high-level pattern:

```cpp
swcme::sep::Interface3D sep(params, spectrum);
auto step = sep.prepare(time_s);       // once per background update

swcme::sep::BackgroundState bg;
auto status = sep.evaluate_background(step, position_m, bg);
```

`BackgroundState` is entirely SI: density `[m^-3]`, proton pressure `[Pa]`,
velocity `[m/s]`, magnetic field `[T]`, magnetic-field magnitude `[T]`,
`div(V)` `[s^-1]`, and Parker focusing length `[m]`.  It also records the exact
prepared model identity and configuration digest.  The 3-D single-point path
delegates to the existing checked batch evaluator with `N=1` and stack scalars,
so an AMPS hot loop does not require a temporary vector or a second status
convention.  After setup and `prepare()` complete, a prepared step and its
owning model may be shared read-only by worker threads.  Each call must still
own distinct destination objects/arrays; model configuration, preparation, and
file output remain setup/coordination operations rather than concurrent
field-query operations. PST04 verifies thread safety, PST07 verifies
direct/adapter equivalence for this complete record, and PST05 verifies that
the state/owner relationship remains defined across moves and destruction.

### Prepared-state immutability (PST01)

A model now has two explicit lifecycle phases.  During initial configuration,
the legacy 1-D fluent setters remain available so existing setup code continues
to compile.  The first **successful**
`prepare_step()` atomically freezes the model configuration.  From that point,
every 1-D setter and 1-D/3-D copy assignment throws
`std::logic_error` before changing any field or model identity.  Failed
validation or numerical preparation does not freeze the model because no valid
state was issued.

The former 1-D `MutableParams()` raw-reference escape hatch is explicitly
deleted.  A reference obtained during setup could otherwise be retained and
used after preparation, bypassing any runtime lock.  Existing code should copy
the read-only `GetParams()` result, edit that standalone value, and pass it to
`SetParams()` before preparation or to `reconfigured()` afterward.

This prevents an old state from changing meaning after construction.  In
particular, changing `sin_theta` can no longer alter the Parker `Bphi` returned
for an already prepared 1-D state.  The same state and fixed query points remain
bitwise identical after every rejected mutation attempt.  The lifecycle flag
is atomic so read-only evaluators may safely inspect a frozen model from the
same threading/MPI contexts already supported by `StepState`; configuration and
preparation remain an exclusive setup phase.

To change a frozen configuration, construct a new owner with
`model.reconfigured(params)`.  The replacement has a fresh `ModelIdentity`, is
independent of all existing states, and remains configurable until its own
first successful preparation:

```cpp
swcme1d::Model model(params);
auto original_step = model.prepare_step(t);  // freezes model

swcme1d::Params changed = model.GetParams();
changed.sin_theta = 0.0;
auto replacement = model.reconfigured(changed);
auto replacement_step = replacement.prepare_step(t);

// model + original_step remain valid and unchanged;
// replacement + replacement_step represent the new configuration.
```

`configuration_locked()` reports the lifecycle state without exposing a way to
reset it.  Copy construction also creates a fresh, initially configurable
owner; it never transfers permission to consume the source model's states.

### Prepared-state model ownership (PST02)

Every `swcme1d::Model` and `swcme3d::Model` now receives a process-unique,
nonzero `ModelIdentity`.  `prepare_step()` stamps that identity into
`StepState::owner_model_identity`.  A state may be shared freely among threads
that evaluate the **same** model instance, but it may not be passed to another
model instance—even if the two instances were constructed with numerically
identical parameters.  This prevents cached geometry, kinematics, Parker, and
shock quantities from one model from being combined with runtime parameters
owned by another model.

Checked APIs reject a foreign or default-constructed state with
`StatusCode::StateModelMismatch` (`STATE_MODEL_MISMATCH`) before validating
other arguments, evaluating physics, clearing destination objects, allocating
mesh/source storage, or opening an output file.  `ModelStatus` then contains:

- `expected_model_identity`: identity of the receiving model;
- `supplied_model_identity`: identity recorded by the supplied state; and
- `has_model_identities=true`, which makes both values part of `summary()`.

The source-compatible value/bool/void geometry and evaluator wrappers cannot
return `ModelStatus`; they call the same ownership guard and throw
`std::runtime_error` on mismatch.  New integrations should use checked paths,
including `shock_acceleration_state_checked()` and the 1-D
`write_tecplot_radial_profile_checked()`, whenever status propagation is
available.  The SEP adapters validate ownership at their outer boundary so
`BackgroundState`, `SEPSourceState`, `SourceSurface`, and optional connectivity
outputs remain unchanged on rejection.

Copy construction creates a new model identity.  Copy assignment can replace
an unprepared model and rotates its identity, but PST01 rejects assignment once
the receiver has successfully prepared a state.  PST03, described below, adds
an independent configuration-snapshot check for foreign or corrupted records.

Correct and incorrect usage therefore look like:

```cpp
swcme3d::Model model_a(params);
swcme3d::Model model_b(params);       // equal Params, different owner
auto step = model_a.prepare_step(t);

auto ok = model_a.evaluate_cartesian_fast_checked(
    step, x, y, z, n, vx, vy, vz, count);       // OK
auto rejected = model_b.evaluate_cartesian_fast_checked(
    step, x, y, z, n, vx, vy, vz, count);       // STATE_MODEL_MISMATCH
```

### Prepared-state configuration ownership (PST03)

Every 1-D and 3-D `StepState` also records a deterministic 64-bit
`configuration_digest`.  The digest covers every public `Params` field in
declaration order, including inactive and deprecated compatibility fields,
data-driven table lengths and values, plus resolved conventions that are not
runtime parameters: configuration-schema version, frame, Parker normalization,
Parker radial polarity, and the proton-only thermal-pressure closure.  Explicit
default values therefore hash identically to implicit defaults, while any
physics-relevant configuration difference changes the snapshot.

The digest uses fixed-width byte serialization and a versioned schema tag; it
does not hash C++ object memory, addresses, padding, locale-formatted text, or a
model identity.  Signed zero and NaN payloads are canonicalized for stable
diagnostics.  It is an allocation-free consistency fingerprint, not a
cryptographic authenticator.

State consumers validate provenance before arguments, allocations, output
mutation, or file opening, using this precedence:

1. A different `owner_model_identity` returns `STATE_MODEL_MISMATCH`.  The
   status includes both model identities and both configuration digests, so
   logs distinguish equal-model misuse from cross-configuration misuse.
2. A matching owner with a different digest returns
   `STATE_CONFIGURATION_MISMATCH`.  PST01 prevents supported model APIs from
   creating this condition; the check remains defense-in-depth for corrupted,
   incompatible, or externally deserialized state records.
3. Only matching owner and configuration snapshots proceed to physics.

For digest diagnostics, `expected_configuration_digest` is the receiving
model's current configuration and `supplied_configuration_digest` is the
snapshot carried by the state; `has_configuration_digests` indicates their
presence.  `ModelStatus::summary()` prints both as fixed-width hexadecimal.
The 3-D model remains construction-time immutable, but it participates in the
same complete digest contract and reports configuration differences on foreign
states.  Use the PST01 replacement path when either dimensional model needs a
different configuration.

### Prepared-state record integrity (PST06)

PST01 prevents model-parameter changes and PST02/PST03 validate provenance,
but those checks alone cannot detect a caller overwriting a public cached value
inside `StepState`.  PST06 therefore gives every successfully prepared 1-D and
3-D record a private 64-bit integrity seal.  The seal covers the owner and
configuration digests, canonical common cache, region/acceleration records,
geometry and Parker caches, Rankine-Hugoniot state, and every retained public
compatibility mirror.

Serialization is explicit and field-by-field.  It never hashes struct memory,
padding, addresses, locale-formatted text, or the seal itself.  Consequently,
normal copy/move construction and assignment preserve an exact valid record
independently of compiler ABI padding.  `integrity_digest()` returns the seal by
value for diagnostics; callers cannot obtain a reference or rewrite it.

Some `StepState` cache members remain publicly visible for source compatibility
with existing diagnostics and demos.  They must be treated as read-only.  If
compatibility code modifies any canonical field or mirror, the next consuming
API recomputes the record digest and returns
`StatusCode::StalePreparedState` (`STALE_PREPARED_STATE`) before evaluating
physics, modifying an output object, allocating result storage, or opening an
output file.  Legacy APIs without a `ModelStatus` return throw
`std::runtime_error` carrying the same status name.

The stale-state diagnostic contains:

- `expected_state_integrity`: the private seal created by `prepare_step()`;
- `computed_state_integrity`: the digest recomputed from the supplied record;
- `has_state_integrity=true`.

Ownership and configuration checks retain precedence.  Corrupting an owner ID
reports `STATE_MODEL_MISMATCH`, while corrupting the configuration tag reports
`STATE_CONFIGURATION_MISMATCH`; all other record corruption reports
`STALE_PREPARED_STATE`.  Do not repair or reseal a rejected state—discard it and
obtain a new state from its owning model.

### Concurrent prepared-state evaluation (PST04)

PST04 turns the documented read-only threading promise into an executable
regression gate.  One 1-D state and one 3-D state are prepared before any
worker starts.  The same owning model and prepared record are then shared by
1, 2, 4, and 8 threads; configuration and the prepared records remain const,
and each invocation writes only to thread-local output storage.

The workload runs scalar AMPS background queries, direct full-field batches,
1-D/3-D source conversion, a 3-D directional shock calculation, and 3-D
Parker-line connectivity.  Separate invalid-input calls execute alongside the
successful work so a hidden shared status buffer would be detected through a
changed code, context, sample index, offending value, or partially written
output.  Connectivity serialization includes all roots and their complete
shock records, covering the allocating update-cadence path as well as the
allocation-free scalar hot path.

For each operation, PST04 explicitly serializes every public status and
numerical field without hashing object padding.  Eight repetitions are run
under forward-interleaved, reverse-interleaved, and operation-grouped schedules
for every required thread count.  Every word must equal the serial oracle
bit-for-bit; there is currently no declared roundoff relaxation.  Both private
prepared-state seals must also remain valid after the stress run.

This guarantee is deliberately scoped.  It does not make model setters,
`prepare_step()`, caller-shared output buffers, or Tecplot writers concurrent.
Those remain externally synchronized setup/output operations.  The ordinary
gate is:

```sh
cd test
make -j
./output/test_swcme --test PST04
```

On a ThreadSanitizer-capable compiler/runtime, instrument the complete test
executable and stop at the first reported race:

```sh
make clean
make CXXFLAGS="-O1 -g -std=c++17 -Wall -Wextra -Wpedantic -fsanitize=thread -fno-omit-frame-pointer"
TSAN_OPTIONS=halt_on_error=1 ./output/test_swcme --test PST04
```

ThreadSanitizer startup failures caused by unsupported container/kernel address
layouts must be recorded as unavailable, not interpreted as a race-free pass.

### Checked output failure propagation (OUT02)

All production 1-D and 3-D Tecplot writers now check the complete buffered
output lifecycle: open, every formatted write, flush, the persistent stream
error indicator, and close.  They share `swcme_output.hpp`, whose
`CheckedTextFile` formats each record into owned memory and then performs an
exact byte-counted write.  This is necessary because a buffered `stdio` write
can appear successful even when the destination rejects the data later during
flush or close, as `/dev/full` does on POSIX systems.

Checked writer APIs distinguish three output failure classes:

- `StatusCode::FileOpenFailure` means that no output handle was acquired;
- `StatusCode::FileWriteFailure` means that a handle was opened but a write,
  flush, stream-error check, or close failed.
- `StatusCode::FileCommitFailure` means that the complete temporary product
  was closed but could not atomically replace the requested destination.

A write failure records `has_io_byte_offset=true` and `io_byte_offset`, the
number of bytes accepted before the first failed operation.  Its stable
`context` identifies the title, variables, zone, data-block/row, flush,
stream-error, or close phase.  Data rows also use `sample_index` for their row
or item number.  Cleanup always attempts close, but the diagnostic is
first-error-wins: a later close failure cannot replace the more useful partial
write location.

The affected status-returning entry points are:

- `swcme1d::Model::write_tecplot_radial_profile_checked()`;
- `swcme1d::Model::write_tecplot_shock_vs_time_checked()`;
- `swcme3d::Model::write_shock_surface_center_metrics_tecplot_checked()`;
- `swcme3d::Model::write_tecplot_dataset_bundle_checked()`; and
- `swcme3d::Model::write_box_face_minX_tecplot_structured_checked()`.

Their source-compatible boolean wrappers now delegate to the same checked
implementation and return `false` for either open or post-open output failure.
New integrations should prefer the checked forms so logs retain the precise
failure phase and byte offset.

The final optional `FileOperations*` parameter on checked writers is intended
for deterministic validation.  Omitting it selects the immutable production
`stdio` backend.  A test backend can fail at an exact accepted-byte count or at
flush/error/close without changing production physics or relying on a specific
filesystem.  OUT02 combines that injected fault matrix with real `/dev/full`
probing at the shared direct-stream layer beneath all five model products.

OUT02 detects and propagates every incomplete write.  OUT03, described next,
adds destination preservation and atomic publication around that checked
stream lifecycle.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT02
```

### Transactional output commit (OUT03)

Every production Tecplot writer now opens an exclusively created staging file
whose name is formed by appending `.swcme-tmp-N` to the destination path.  The
staging file is therefore in the destination directory and on the same
filesystem.  The writer sends all records through the OUT02 checks, flushes,
examines the stream error state, and closes the staging handle before calling
`rename()` to install it.  On POSIX systems—the supported SWCME production
environment—the same-directory rename atomically replaces an existing regular
file, so readers observe either the previous complete product or the new
complete product, never a partially written zone.

The commit sequence is:

1. validate model state and all caller inputs before filesystem access;
2. exclusively create a unique same-directory staging file;
3. write and validate the complete product through OUT02;
4. close the staging handle successfully;
5. atomically rename the staging file over an absent or regular destination;
6. remove the staging file after any write, close, or commit failure.

Commit is deliberately refused when the destination is a directory, symbolic
link, device, FIFO, or other non-regular object.  This prevents a privileged
validation run from replacing objects such as `/dev/full`.  If staging cannot
be created, checked writers return `FILE_OPEN_FAILURE`; stream failures retain
`FILE_WRITE_FAILURE`; and a failed final rename returns
`FILE_COMMIT_FAILURE`.  The latter includes `io_byte_offset`, which in this
case is the complete staged-product size.  A cleanup failure never masks the
earlier write or commit diagnostic, although an operating-system refusal to
remove may leave the private staging file for diagnosis.

The optional `FileOperations` validation table now includes `commit` and
`remove` callbacks.  Production callers continue to omit the table.  Legacy
boolean writer APIs use the same transaction and return `false` for any open,
write, close, or commit failure.

OUT03 provides atomic namespace publication, not power-loss durability: it
does not currently call `fsync()` on the file and containing directory.  Runs
that require recovery across a host crash should add that separate durability
policy at the campaign-storage layer.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT03
```

### Model-domain output preflight (OUT05)

All five Tecplot products now complete their model-domain validation before
selecting an output backend or creating an OUT03 staging file.  This closes the
remaining gap between checked evaluation and checked publication: atomic commit
kept an existing destination safe, but a bad late sample could still perform
filesystem work and, in the shock-history writer, prepare only part of the
requested series before failing.

The 1-D radial-profile writer scans every radius and supplied field value.
Non-finite radii return `NONFINITE_INPUT`, radii below
`solarwind::MIN_RADIUS_M` return `OUTSIDE_MODEL_DOMAIN`, and non-finite
precomputed fields return `NONFINITE_RESULT`; each status carries the first
bad row in `sample_index`.  The shock-history writer prepares and validates the
complete immutable time series in memory before opening output, so an invalid
late data-driven time also identifies its requested time and row without an
output side effect.

For every structurally valid box (see OUT04 below), the 3-D writers use one
shared coordinate generator for preflight and emission.
Surface nodes are checked explicitly, the bundle scans all structured volume
points plus its min-X face, and the standalone face scans the exact grid it will
write.  Finite points below the radius floor are rejected as
`OUTSIDE_MODEL_DOMAIN`.  The boundary itself is inclusive.

OUT05 preserves the established ordering of safeguards:

1. reject a foreign, reconfigured, or corrupted prepared state;
2. reject null and malformed non-box arguments;
3. validate surface mesh/metrics through OUT06 where present;
4. validate the complete BoxSpec through OUT04 where present;
5. preflight every requested output sample and generated coordinate;
6. create a private staging file and run the OUT02 checked-write lifecycle;
7. publish the complete product through the OUT03 atomic commit.

Legacy boolean writers delegate to their checked companions, so they return
`false` under the same no-open domain preflight.  Checked callers additionally
receive the precise status, offending value, context, and sample index.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT05
```

### Box specification validation (OUT04)

`BoxSpec` now has one enforced structural contract shared by the bundle,
standalone min-X face, and `default_apex_box()` factory:

- `cx`, `cy`, `cz`, `hx`, `hy`, and `hz` must be finite;
- half extents must be nonnegative (zero remains valid for collapsed diagnostic
  grids);
- `Ni`, `Nj`, and `Nk` must each be at least two, including `Ni` for the
  standalone face so a `BoxSpec` is valid independent of its consumer;
- each center-minus-extent bound, center-plus-extent bound, and doubled span
  must be representable as a finite `double`; and
- the complete `Ni*Nj*Nk` cardinality must be representable as `size_t`, using
  checked multiplication rather than a potentially wrapped product.

Non-finite members or derived bounds/spans return `NONFINITE_INPUT` with the
offending value.  Negative extents, undersized dimensions, and cardinality
overflow return `INVALID_CONFIGURATION`.  Validation is performed after the
prepared-state ownership/integrity gate but before OUT05 point traversal,
physics evaluation, output-backend selection, or staging-file creation.
Legacy boolean writers delegate to the checked implementations and return
`false` for the same invalid specifications.  The default-box factory throws
the corresponding status exception rather than returning an unusable box.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT04
```

### Mesh output validation (OUT06)

The surface and four-zone bundle writers now treat `ShockMesh` and
`TriMetrics` as one validated output record rather than unrelated arrays that
only need matching lengths.  Before model-domain traversal or output access,
OUT06 enforces:

- at least three vertices and one triangle;
- exact `Nv` lengths for every nodal array and exact `Ne` lengths for all three
  connectivity arrays;
- finite nodal positions, normals, compression, and normal shock speed;
- unit nodal normals within `1e-10`, compression `rc >= 1`, and
  `Vsh_n >= 0`;
- one-based, in-range, distinct indices for every triangle;
- the existing scale-relative nondegeneracy and outward-winding checks from
  `compute_triangle_metrics()`; and
- one unambiguous metric mode: every `TriMetrics` vector empty to request
  canonical computation, or every vector complete and numerically consistent
  with a fresh canonical derivation from the supplied mesh.

Complete caller-supplied metrics are compared field-by-field against the
canonical area, normal, centroid, mean compression, and mean normal speed.
This detects stale finite metrics after coordinates, connectivity, `rc`, or
`Vsh_n` change—defects that size and finiteness checks cannot find.  The scaled
tolerance is `1024*epsilon`, substantially tighter than the nine-digit Tecplot
serialization.  Partial metrics are rejected rather than silently discarded.

Structural, topology, physical-invariant, degenerate-cell, partial, and stale
records return `INVALID_MESH`; non-finite model fields return
`NONFINITE_RESULT`.  Where a node or connectivity/metric row is known, its
index is carried in `sample_index`.  Both checked writers reject before
selecting `FileOperations` or creating an OUT03 staging file.  Legacy writers
delegate to the checked path and return `false`.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT06
```

### Independent output parsing (OUT01)

The three production 3-D Tecplot products now use one explicit 25-variable
schema whose names include their physical units.  The exact ordered tokens are
`X[m]`, `Y[m]`, `Z[m]`, `n[m^-3]`, the three velocity components in `[m/s]`,
the three magnetic-field components in `[T]`, `divVsw[s^-1]`, the
dimensionless compression and normal fields marked `[-]`, shock speeds in
`[m/s]`, triangle area in `[m^2]`, and centroids in `[m]`.  This replaces the
former unitless 3-D header names; downstream scripts should select the new
unit-qualified names.

OUT01 writes a small deterministic surface-only file, four-zone dataset
bundle, and standalone min-X face through the normal transactional filesystem
backend.  A validation-only parser then reads the committed bytes using its own
Tecplot grammar and a separately declared expected schema.  It does not import
the writer's variable literal or output helpers, so a production-side format
change cannot silently change the reference side of the test.

The parser verifies quoted titles, variable names, units and order; exact zone
kinds and declarations; BLOCK variable locations and per-block counts; POINT
row widths; structured dimensions; node and element counts; one-based,
in-range, distinct triangle indices; numerical finiteness; and exact
end-of-file consumption.  Parsed connectivity is compared with the requested
mesh, and representative surface, volume, and face values are compared with
direct model evaluations within the ten-significant-digit `%.9e` formatting
precision.  The standalone face must also reproduce the face zone embedded in
the bundle.

Sensitivity probes corrupt one contract at a time—an added trailing record, a
changed variable unit, an out-of-range triangle index, and an extra POINT
column—and require the independent parser to reject each mutation.  OUT01 is a
post-commit consumer validation; the earlier OUT06, OUT04, OUT05, OUT02, and
OUT03 gates continue to own preflight, write detection, and publication safety.

Run the gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT01
```

### Demonstration program execution (OUT07)

The three user-facing examples are now first-class, warning-checked build
products.  From `test/`, `make demos` builds `output/demo1d`,
`output/demo3d_1`, and `output/demo3d_2` with the same compiler flags and 3-D
production object used by the validation executable.  `make demo-run` builds
the complete gate and runs OUT07.

The two 3-D examples no longer construct visualization boxes whose faces pass
through the solar origin.  They use `default_apex_box()` to create a validated
12-by-12-by-12 box centered near the shock apex, so every volume and min-X-face
sample is inside the supported radius.  Their demonstration meshes use 24
polar intervals and 48 periodic azimuthal nodes: enough to exercise the unique
apex, finite SSE boundary, connectivity, and area metrics without the former
multi-minute, very large output.  The bundle name is now
`sse_apex_bundle_tecplot.dat`, and comments describe its actual four zones in
production order.

Output failures are no longer advisory.  Both 3-D examples use the checked
bundle API and terminate nonzero with the full status diagnostic if writing
fails.  CSV streams check open and delayed close state.  The auxiliary point
clouds in `demo3d_2` now use unit-qualified variables and the shared checked,
transactional text-output layer.  The extended demo's strength table also
reports the production ideal-MHD shock result—compression, fast Mach number,
shock-frame speed, `theta_Bn`, magnetic amplification/rotation, and normalized
conservation residuals—instead of inferring Mach number from the obsolete
gas-dynamic compression proxy.

OUT07 executes each binary in a separate freshly emptied directory with a
fixed C locale and captured stdout/stderr.  It requires exit code zero, empty
stderr, the exact declared artifact manifest, and no transaction staging
files.  Every Tecplot product is parsed through the independent OUT01 grammar;
every CSV is checked for its exact header, field count, finite values, final
newline, 865-row five-minute time axis, and 72-hour endpoint.  Surface and box
sizes are independently derived from the documented demo resolution.  Passing
temporary products are removed; failed runs remain under
`test/output/OUT07_demo_runs_<run-id>` with their logs for diagnosis.  The
process-and-start-time identifier keeps simultaneous validation campaigns
isolated from one another, including in containers that reuse process IDs.

Run the complete demonstration gate with:

```sh
cd test
make demo-run
```

or run it after an existing build:

```sh
./output/test_swcme --test OUT07
```

### AMPS adapter equivalence (PST07)

PST07 treats `swcme_sep_interface.hpp` as a transparent integration boundary,
not a second model.  For representative 1-D and 3-D prepared states it compares
the adapter with direct checked SWCME calls at identical times, positions, and
SI units.  The comparison covers success and domain-failure status, prepared
model/configuration identity, density, proton pressure, velocity, Parker-field
components and magnitude, `div(V)`, local magnetic focusing length, shock
position/normal, compression, `theta_Bn`, fast Mach number, shock speed,
upstream state, source weighting, and a five-node relativistic energy spectrum.

Parker path length and focusing now live beside the production Parker field in
`swcme_solarwind.hpp`.  The 3-D connectivity solver delegates its arc-length
calculation to that common helper, while both AMPS adapters use the common
pressure and focusing helpers.  Thus the adapter assembles records and checks
status but contains no separate pressure, path, focusing, shock, or spectrum
physics.  `BackgroundState` records the exact prepared model identity and
configuration digest.  Observer-connected `SEPSourceState` records the selected
Parker path length; all source records include upstream pressure and local
focusing, with unavailable non-observer path lengths represented by `NaN` and
serialized as `NA` rather than an ambiguous zero.
The resolved manifest records this expanded layout as SEP source contract
version 2 so older campaign readers cannot silently interpret shifted columns.

Run the dedicated equivalence gate with:

```sh
cd test
make -j
./output/test_swcme --test PST07
```

PST07 follows OUT07 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

### Prepared-state lifetime contract (PST05)

`StepState` is a self-contained value record. Neither the 1-D nor 3-D record
contains a pointer or reference to its preparing `Model`; all Parker, density,
kinematic, geometry, region, acceleration, and shock data needed by an
evaluation are owned by value. A state may consequently be copied, moved,
archived in memory, inspected, and destroyed after its original model has gone
out of scope without accessing released storage.

Self-contained storage does not make a state an independent evaluator. Every
checked physics or AMPS-adapter call still requires a live logical model owner.
Constructing a new model with numerically identical parameters creates a new
`ModelIdentity`, so it rejects a state whose original owner was destroyed with
`STATE_MODEL_MISMATCH` before modifying caller-owned output. This is the
documented use-after-owner behavior; SWCME never attempts to find or silently
rebuild the destroyed owner.

Move construction and move assignment are different from copying: they
relocate the same logical owner. The model identity and its prepared/frozen
phase transfer to the destination, and states prepared before the move remain
valid there with bitwise-identical results. The moved-from object receives a
fresh identity and is unlocked; it is safe to destroy or assign and cannot
consume the transferred states. Its remaining parameter value is the normal
C++ valid-but-unspecified result of a move and must not be used as a configured
physics model unless it is assigned a complete configuration again.

Move assignment into a destination that has already prepared a state remains
forbidden by PST01. It throws `std::logic_error` before either source or
destination changes, preserving both owners and all their existing states.
Model move construction is conditionally `noexcept` from `Params`, and the
current value-only parameter bundles satisfy that condition. Standard
containers therefore relocate models and adapters through the ownership-
transferring move path instead of the independent-owner copy path.

The same rules apply to `Interface1D` and `Interface3D`: adapter copies own new
models and reject source states, whereas adapter moves transfer the embedded
model identity. Asynchronous evaluation is supported only while the relocated
owner remains alive, with the model and prepared state shared read-only and
each task owning separate output storage. No contract permits a worker to
retain a model reference beyond the model's lifetime.

Run the ordinary and AddressSanitizer lifetime gates with:

```sh
cd test
make -j
./output/test_swcme --test PST05
make pst05-sanitize
```

`pst05-sanitize` compiles a separate fully instrumented executable with
`-fsanitize=address` and fail-fast behavior, then runs PST05. It does not reuse
or alter normal incremental objects. LeakSanitizer is disabled by default
because it cannot run under ptrace-based CI/container supervisors. On an
untraced host, enable it with
`make PST05_ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 pst05-sanitize`.
PST05 follows PST07 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through `@ALL`.

### State-ownership performance (PST08)

PST08 makes prepared-state safety a bounded-cost contract. Public APIs still
validate model identity, the complete configuration digest, and the PST06
state-integrity seal before inspecting arguments or modifying outputs. The
implementation now guarantees that this complete authentication occurs once
per outer public call—never once per nested evaluator, divergence stencil, or
batch sample.

Two changes enforce that rule. First, each model caches its complete
configuration digest during construction or guarded setup. The 1-D fluent
setters refresh the cache before preparation; the first successful preparation
then freezes both parameters and digest. Copy and move operations preserve the
PST02/PST05 identity rules while maintaining the matching cached digest.
Validation therefore compares one scalar instead of serializing the complete
parameter bundle in every particle query.

Second, composite operations use private `*_after_validation` kernels. These
kernels are inaccessible to callers and are entered only after an outer public
method has authenticated the state. The 1-D full-field evaluator no longer
re-enters the checked fast evaluator, the 3-D field loop no longer revalidates
inside directional shock/geometry calls, and Cartesian divergence stencils no
longer hash the same record for each velocity sample. The AMPS adapters build a
candidate output locally and delegate authentication to one checked Model call;
state failures still leave the caller's prior record unchanged.

The dedicated benchmark separates `prepare_step()` from evaluation and reports
median and nearest-rank p95 timing after warm-up. It measures complete state
validation, direct and AMPS scalar queries, and representative 1-D/3-D batches.
An emulated pre-PST08 topology repeats the exact current validation at the
locations removed by the correction, providing an on-host before/after cost
comparison without exposing an unsafe unvalidated API. A test-only global
allocation probe covers scalar, array, sized, and aligned C++17 allocation and
requires zero allocations in warmed validation and direct/AMPS hot paths.

The declared acceptance budgets are intentionally conservative for shared CI:

- complete validation p95 is at most 10 microseconds;
- 1-D direct and AMPS scalar p95 are at most 20 and 25 microseconds;
- 3-D direct and AMPS scalar p95 are at most 500 and 600 microseconds;
- a 16,384-point 1-D batch completes within 5 milliseconds p95;
- a 256-point 3-D batch completes within 100 milliseconds p95;
- one validation is at most 5% of either representative batch runtime;
- corrected 1-D direct and AMPS medians improve materially over their emulated
  duplicate-validation topology; and
- warmed ownership, direct, and AMPS hot paths allocate no memory.

Run PST08 through the ordinary suite or its pinned optimized target:

```sh
cd test
make -j
./output/test_swcme --test PST08
make pst08-performance
```

The pinned target compiles only the production 3-D implementation and PST08
with C++17, `-O3`, `-DNDEBUG`, warnings, and pthread support, then prints the
complete timing record. Absolute nanosecond values are hardware- and load-
dependent; the budgets and same-process before/after ratios are the release
gate. PST08 follows PST05 in `SMOKE`; other profiles include it through `@ALL`.

### Strict-warning writer build (OUT08)

OUT08 turns compiler diagnostics for the output boundary into a release gate.
The dedicated target builds the shared status/output headers, compiled 3-D
writers, header-only 1-D writer, and all three demonstrations twice: once with
`-O0 -g3` and once with `-O3 -DNDEBUG`. Both variants enable `-Wall`,
`-Wextra`, `-Wpedantic`, `-Wformat=2`, `-Wformat-security`, `-Wconversion`,
`-Wsign-conversion`, and `-Wshadow`; `-Werror` makes every selected diagnostic
fatal, with `-Werror=format-security` retained explicitly in the command line
so the security requirement remains visible in CI logs.

The checked variadic formatter now carries a GCC/Clang `printf` format
annotation. Consequently the compiler verifies every production format literal
against its argument types instead of treating `CheckedTextFile::print()` as an
opaque project function. Records without substitutions can use the templated
`write_literal()` path, which accepts a compile-time-sized character array and
never interprets percent tokens. Both paths share one checked raw-record helper,
so partial-write byte accounting and first-error-wins behavior remain identical.

The strict optimized build exposed possible use of uninitialized scalars in the
3-D volume and face writers. Those writers now initialize every output field,
call status-returning field and divergence evaluators, remap a failure to the
global output-row index, and cancel the staging transaction before returning.
The fix therefore addresses the underlying control-flow ambiguity rather than
suppressing the warning. Signed mesh connectivity is likewise validated before
one explicit conversion to `std::size_t`; array-axis loops and demonstration
time/grid conversions now use explicit, type-correct boundaries. OUT08 contains
no warning suppression.

Run the registered record-semantics assertion and strict compiler gate with:

```sh
cd test
make -j
./output/test_swcme --test OUT08
make out08-strict CXX=g++
# Run the same compiler-independent target in a Clang CI job when available:
make out08-strict CXX=clang++
```

`out08-strict` creates isolated `output/out08-debug/` and
`output/out08-optimized/` products, compiles every target from source rather
than reusing normal objects, and runs the same literal/format assertion in both
configurations. The runtime assertion verifies that `%` sequences remain
verbatim through `write_literal()` and that a correctly typed formatted record
still follows the checked open/write/flush/error/close lifecycle. Acceptance
requires both builds and both assertions to pass without a diagnostic. OUT08
follows PST08 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include its registered
runtime half through `@ALL`, while CI/release jobs must additionally invoke the
compile-time `out08-strict` target.

### `SEPSourceState`

`SEPSourceState` is the stable transport-facing source record.  It contains:

- explicit `active`, `connection_evaluated`, and `connected` flags;
- launch-relative time and deterministic source ID;
- source position `[m]`, outward normal, patch area `[m^2]`, active surface
  area, area fraction, and relative patch weight;
- compression ratio, `theta_Bn`, fast Mach number and normal shock speed;
- upstream density, pressure, and magnetic-field magnitude;
- local Parker focusing length and, for observer cobpoints, field-line path
  length;
- the DSA phase-space slope `q` and derived intensity indices;
- the complete source-spectrum configuration and normalization convention.

Surface sources are generated with
`Interface3D::build_shock_surface_source()`.  One record is returned per
validated triangular shock-mesh cell.  Only cells containing a physical fast
shock in `SOURCE` mode receive nonzero source weight.  Their area fractions are
normalized over the physical active area, so a spatially uniform relative
source is sampled according to actual surface area rather than mesh index.

An observer-connected source is obtained with
`Interface3D::source_at_observer_cobpoint()`.  It uses the production Parker
field-line/shock intersection and routes the selected cobpoint direction
through the same canonical source-state path used by shock-surface patches.
`NO_CONNECTION` is an expected status, distinct from a numerical or shock-solver
failure.

For the exact spherical/+X reduction, `Interface1D::source_at_shock()` and the
3-D direction adapter emit byte-identical deterministic source records.  This is
protected by `SEP03` in addition to the lower-level `1D3D03` acceleration test.

### Spectrum and normalization convention

In `SOURCE` mode the shock module supplies

```text
f(p) proportional to p^(-q),      q = 3 r_c/(r_c-1).
```

For an isotropic distribution the SEP adapter uses the exact relativistic
intensity shape

```text
J(E) / J(E_ref) = [p(E) / p(E_ref)]^(2-q),
p c = sqrt[K (K + 2 m c^2)].
```

The public energy inputs are MeV; momentum is calculated in SI.  Helpers are
provided for rigidity in GV and for converting differential intensity between
`(cm^2 s sr MeV)^-1` and `(m^2 s sr J)^-1`.

Two normalization modes are explicit:

- `RELATIVE_ONLY` (default) supplies a dimensionless source shape/patch weight
  and deliberately does not claim an absolute injection flux;
- `REFERENCE_DIFFERENTIAL_INTENSITY` requires a positive physical
  `J(E_ref)` in SI and produces a dimensional differential intensity.

Thus a relative validation source can never be silently exported with physical
flux units.  In `RESOLVED_COMPRESSION` mode the prescribed source is inactive;
requesting its spectrum returns `SOURCE_INACTIVE` rather than a fabricated DSA
spectrum.

`swcme::sep::source_csv_header()` and `serialize_source_csv()` provide a stable,
audit-oriented CSV representation.  They are intended for regression and
AMPS/SWCME handoff checks, not as a mission archive standard.

### Reference consumer

`tools/sep_reference.cpp` is a small executable that consumes only the public
SEP interface.  It can produce a 1-AU cobpoint/source time history, print the
resolved configuration manifest, or emit machine-readable `KEY=value` probe
records for campaign sweeps:

```sh
make -C test tools
./test/output/sep_reference --print-manifest
./test/output/sep_reference --start-hours 0 --end-hours 72 --step-hours 1 \
    --output reference.csv
./test/output/sep_reference --probe-time-hours 24 --v0-kms 1500
```

It is intentionally a reference consumer rather than an alternative physics
implementation.

## Automated validation campaigns

`test/run_tests.py` is the higher-level validation/campaign manager.  It always
controls the single C++ `test/output/test_swcme` executable; Python does not
reimplement the deterministic physics tests.

```sh
cd test
python3 run_tests.py --list
python3 run_tests.py --profile SMOKE
python3 run_tests.py --profile ROUTINE
python3 run_tests.py --all                  # FULL
python3 run_tests.py --test SEP03
python3 run_tests.py --profile EVENT --event-config event_config.example.json
```

Profiles are stored in `test/profiles/`:

- `SMOKE` is a short development gate covering prepared-state safety, checked
  output failure propagation, transactional commit, model-domain preflight,
  BoxSpec and mesh-output validation, independent parsing, executable
  demonstrations, state-ownership performance, strict writer record semantics,
  configuration, core shock, connectivity, divergence, and SEP-interface
  integration;
- `ROUTINE` runs the broad deterministic suite while excluding the slowest
  stochastic/multi-root stress cases;
- `FULL` runs the complete registered C++ suite and exports the default SEP
  reference history unless disabled;
- `EVENT` begins with the FULL verification gate and then executes the
  event-specific JSON analysis specification.

Every campaign writes `manifest.json`, `summary.json`, `summary.csv`, and one
log per C++ test.  The manifest records the selected tests and seed, git
commit/branch/dirty state, compiler path/version/flags, host/Python information,
visible MPI and OpenMP environment, a source-tree SHA-256, the complete resolved
SWCME+SEP configuration and its SHA-256, and the event JSON/hash when present.
This prevents an event comparison from being detached from the exact model
configuration that produced it.

EVENT JSON supports Cartesian parameter sweeps with placeholder substitution,
regex metric extraction and bounds; log-log convergence-order fits; keyed or
row-aligned CSV comparisons with absolute/relative tolerances; and optional
Matplotlib plots.  Commands may use `{seed}`, `{root}`, `{test_dir}`, and
`{output_dir}` in addition to sweep parameter names.  See
`test/event_config.example.json` for an executable example and `test/README.md`
for the full schema.

Campaign exit codes are deterministic: `0` = all required work passed, `1` = a
validation/reference/event analysis failed, `2` = command/configuration error,
and `3` = build failure.  `test/python/test_run_tests.py` regression-tests the
profile expansion, convergence fitting, sweep engine, CSV comparison/event
aggregation, manifest/report production, and command-error exit code.

## Current 3-D shock geometries

`swcme3d::ShockShape` currently provides:

- `Sphere` — a Sun-centered spherical verification geometry;
- `Ellipsoid` — a Sun-centered ellipsoid whose three axes scale self-similarly
  with the shock apex distance; and
- `SSE` — the recommended **and default** finite self-similar-expansion spherical cap.

`ConeSSE` is retained as an enum alias for source compatibility but now has the
same corrected semantics as `SSE`.  It no longer means the old
`R=R_apex*cos(theta)^m` cosine-cap model.  The legacy `flank_slowdown_m`
parameter is retained in `Params` only so existing source/input code continues
to compile; it is ignored for the corrected SSE geometry.

### Finite SSE geometry

For apex distance `R_apex` and half width `lambda`, the generating sphere is

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda),
```

where `c` is the heliocentric distance of the sphere center along the CME axis
and `a` is its radius.  A ray separated by `alpha` from the CME axis intersects
the outward cap at

```text
R(alpha) = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2)
```

for `alpha <= lambda`.  The ray is tangent at `alpha=lambda`; no surface exists
for larger angular separation.  The exact outward normal is

```text
n_hat = (R e_r - c e_CME) / a.
```

`Model::shape_radius_normal()` returns a boolean surface-existence flag.  When
it returns `false`, the radius and normal outputs are zero and are not physical
values.  `Model::diagnose_direction()` follows the same convention and returns
`rc=1`, `V_sh,n=0` for an absent finite surface.  The Cartesian evaluators use
this flag internally and return the undisturbed ambient wind/Parker field
outside the SSE angular support.

### Self-similar flank speed

All supported shapes use the same geometric normal-speed relation.  If the
shape scales with apex distance so that `R(u,t)=f(u) R_apex(t)`, then

```text
dR/dt   = V_apex R/R_apex
V_sh,n  = V_apex (R/R_apex) (e_r dot n_hat).
```

This replaces the legacy independent cosine flank-speed factor and corrects the
ellipsoid flank speed, which previously used the full apex speed away from the
apex.


## Shared CME/shock-apex kinematics

The 1-D and 3-D models now obtain the shock-apex radius and speed from the same
`swcme_kinematics.hpp` implementation.  This removes the former divergence in
which the 1-D slow-CME branch clipped `V0-Vsw` to zero while the 3-D branch used
the fast-CME formula for negative speed differences and divided by `Gamma` at
`Gamma=0`.

Three kinematic modes are supported:

- `swcme::kinematics::Mode::Ballistic` uses
  `R=R0+V0*t`, `V=V0` exactly;
- `Mode::DBM` uses the sign-aware constant-background drag-based solution; and
- `Mode::DataDriven` uses monotone PCHIP interpolation of a supplied
  height-time table and returns the derivative of that same interpolant as the
  apex speed.

For DBM, define `DeltaV0=V0-Vsw` and `a=abs(DeltaV0)`.  The common solution is

```text
DeltaV(t) = DeltaV0 / (1 + Gamma a t)
R(t)      = R0 + Vsw t
            + sign(DeltaV0) log(1 + Gamma a t) / Gamma.
```

The absolute value in the drag denominator is essential: fast CMEs decelerate
toward the ambient wind and slow CMEs accelerate toward it.  `Gamma=0` is an
explicit ballistic branch, so there is no division by zero and no physics-level
`finite_or()` fallback.  A short series for `log1p(x)/x` is used at very small
`x=Gamma*a*t` so the nonzero-drag solution remains continuous with the exact
ballistic limit.

### Data-driven kinematics

`Params::data_time_s` contains strictly increasing times in seconds and
`Params::data_radius_Rs` contains nondecreasing apex radii in nominal solar
radii.  PCHIP was selected because it passes exactly through the supplied
height-time knots while preserving monotonicity and avoiding cubic overshoot.
The local derivative is used as `V_sh`, so the reported speed is kinematically
consistent with the radius curve.

The default extrapolation policy is
`swcme::kinematics::ExtrapolationPolicy::OutsideTime`: a query before the first
or after the last knot is rejected explicitly.  `Ballistic` continuation may be
selected when an explicit endpoint continuation assumption is desired; it uses
the endpoint PCHIP derivative and does not silently extrapolate the cubic.

Both dimensional wrappers convert their public parameters to the same SI
`swcme::kinematics::Config`.  Therefore identical kinematic inputs are required
to produce identical apex radius and speed to roundoff.  The recommended
default DBM reference radius is now `20 R_s`, reflecting the intended use of
this simple drag model in the drag-dominated heliosphere rather than at
`~1.05 R_s`.  Event-specific calculations may still set another radius when
there is a documented physical justification; observationally constrained
`DataDriven` mode is preferred when height-time measurements are available.

The deterministic kinematics validation block is `KIN01`-`KIN08`; see
`test/README.md` for individual purposes and acceptance criteria.

## Observer-to-shock magnetic connectivity and cobpoint tracking

The 3-D model now provides an analytical Parker-field-line connectivity solver
through `Model::observer_connectivity()`.  The solver is designed for the
controlled upstream-Parker experiment used by the SEP study: it traces the
observer's nominal Parker line inward and intersects that line with the same
production shock geometry used by `shape_radius_normal()` and the same local
shock physics used by `shock_state_direction()`.

For the production Parker field

```text
B_phi/B_r = -Omega r sin(theta) / V_sw,
```

a field-line tangent satisfies

```text
r sin(theta) dphi/dr = B_phi/B_r,
```

so `dphi/dr=-Omega/V_sw`.  The exact observer-anchored field line is therefore
constructed by rotating the observer radial direction about the configured
solar-rotation axis by

```text
Delta phi = -Omega (r-r_obs) / V_sw.
```

`Params::solar_rotation_rate_rad_s` now explicitly carries the rotation rate
used by both the Parker magnetic field and connectivity mapping.  Its default
is the existing SWCME solar-rotation convention.  Setting it to zero provides
the exact radial-field limit used by `CON01`; the normal production default is
unchanged.

`ConnectivityState` retains every geometrical field-line/shock intersection in
increasing radial order.  Each `ConnectivityRoot` contains the Cartesian
cobpoint, shock radius residual, analytical Parker path length to the observer,
and the complete `LocalShockState`.  The default selected cobpoint is the
outermost root, i.e. the first shock surface encountered when tracing inward
from the observer.  Retaining all roots makes this choice explicit and keeps
the infrastructure usable for future non-convex geometries.

The root search is deliberately robust to connection boundaries.  It combines
radial scanning, bisection of sign-changing roots, local minimization of the
surface residual to detect tangent roots that do not change sign, and explicit
refinement of finite-SSE surface-validity transitions.  No artificial SSE
flank is introduced when the Parker line remains outside the configured cap.

The Parker path length returned with a cobpoint is analytical.  With

```text
k = Omega sin(theta) / V_sw,
```

SWCME integrates

```text
ds/dr = sqrt(1 + (k r)^2)
```

in closed form.  The exact zero-rotation/polar limit is the radial distance.
This length, rather than radial separation, is the quantity intended for
field-aligned SEP transport timing.

`Model::observer_connectivity_history()` evaluates a stationary observer at a
requested set of times.  Every time step is solved independently from the
production kinematics, Parker line, geometry, and shock state; the history
contains no hidden hysteresis.  Consequently connection onset/loss and
cobpoint motion can be interpreted as model physics rather than state retained
by the tracker.

The deterministic connectivity validation block is `CON01`-`CON08`; see
`test/README.md` for the individual fixtures and acceptance checks.

## Centralized configuration validation and unit handling

SWCME now separates **unit conversion** from **physical admissibility**.
`swcme_units.hpp` is the single production source for conversions between the
public heliophysics units and SI, including km/s, nT, cm^-3, km^-1, AU, solar
radii, hours, and degrees.  Conversion functions are deliberately pure: for
example, `0 km/s` converts to `0 m/s`; it is not silently replaced by a
positive speed.

`swcme_config.hpp` provides the common validation contract used by both the
1-D and 3-D models.  `Model::validate()` is side-effect free and returns every
invalid field with a structured code and requirement.  `prepare_step()` calls
the same validator before any basis normalization, unit conversion, kinematic
evaluation, or field/shock calculation.  Invalid input therefore fails once at
setup rather than being clipped into a plausible-looking state.

Common checks include positive solar-wind speed, reference density and
temperature; non-negative magnetic field and DBM drag coefficient; `gamma>1`;
`sin(theta)` in `[0,1]`; valid region/smoothing parameters; and well-formed
DATA_DRIVEN tables. Smoothing widths are validated against their adjacent
layers even when inactive in the selected mode, so an accepted stored
configuration contains no dormant invalid value. The 3-D interface additionally
validates non-zero finite CME and solar-rotation axes, positive ellipsoid axis
ratios, non-negative solar rotation rate, and SSE half width in `(0,pi/2]`.

The production models now use the centralized conversion helpers when preparing
their SI state.  This fixes the prior 1-D `V_sw` path that used
`max(1,V_sw*1000)` and therefore changed `0 km/s` into `1 m/s`.  Zero wind speed
is now converted exactly and rejected by the configuration validator because a
positive wind speed is required by the Parker/DBM baseline.

See `CONFIGURATION_UNITS_FIX_NOTES.md` and the `CFG01`-`CFG03` sections in
`test/README.md` for the full contract and validation coverage.

### Smoothing-width policy (CFG03)

CFG03 removes the former silent 90-percent cap from region construction. Public
widths are total, symmetric, self-similar fractions specified as AU at a 1-AU
shock. For sheath fraction `f_s` and ejecta fraction `f_e`, centralized
configuration validation now requires

```text
w_shock <= 0.90 f_s,
w_LE    <= 0.90 min(f_s, f_e),
w_TE    <= 0.90 f_e.
```

Equality is accepted; the next representable floating-point value above any
limit is rejected with `OUT_OF_RANGE` and the exact public field name. Negative
and non-finite widths retain their existing `NEGATIVE` and `NON_FINITE`
classifications. The rule applies to every stored width, including widths
inactive under `SHOCK_ONLY/SOURCE`, ensuring a later setup-phase mode change
cannot activate a value that validation previously allowed only because it was
dormant.

After validation, `swcme::regions::make_boundaries()` applies each accepted
fraction directly:

```text
W_effective(u) = w_requested R_sh(u).
```

There is no `min`, maximum, adjustment status, or hidden fallback. Thus the
effective width in a prepared 1-D state, at the 3-D apex, and at every local
shock-surface direction agrees with the recorded request to floating-point
roundoff. The public `smoothing_fraction_limits()` helper is the single source
for the 90-percent margins used by both validation and tests.

Run the focused gate with:

```sh
cd test
make -j
./output/test_swcme --test CFG03
```

The test covers zero and ordinary widths, simultaneous exact-limit widths,
one-ULP violations for each field, NaN and infinities, multi-field conflicts,
inactive oversized values, exact local 1-D/3-D scaling, non-overlapping
transition intervals, and finite half-blends. CFG03 follows OUT08 in `SMOKE`;
the other profiles include it through `@ALL`.


## Shared common physics core for 1-D and 3-D

SWCME now prepares dimensionality-independent solar-wind and apex-kinematic
physics through two shared production components:

- `swcme_solarwind.hpp` owns the Leblanc density coefficients and normalization,
  Parker radial-field normalization, Parker scalar components, the Cartesian
  Parker vector for an arbitrary solar axis, and the current proton thermal
  pressure closure; and
- `swcme_core.hpp` converts one common public-unit configuration to SI, prepares
  the shared solar-wind cache, and evaluates the shared ballistic/DBM/data-driven
  apex kinematics.

Both `swcme1d::StepState` and `swcme3d::StepState` retain legacy mirror fields
(`C2/C4/C6`, `Br1AU_T`, `k_AU`, `V_sw`, apex radius/speed) so existing source
continues to compile.  Those fields are no longer independently calculated; they
are copied from `StepState::common`.  This makes the common prepared state the
authoritative source while preserving the current public interfaces.

The dimensional wrappers now differ only where geometry genuinely differs.  The
1-D model supplies a fixed ray `sin(theta)` to the common Parker-component
function.  The 3-D model derives local latitude from the solar rotation axis and
radial direction and asks the same common Parker implementation to construct the
Cartesian field.  Both models call the same common Leblanc density evaluator and
the same proton-pressure closure before entering the already shared ideal-MHD
shock solver.

The new `1D3D01` and `1D3D02` validation tests are permanent guards against a
return of duplicated physics.  `1D3D01` compares the canonical common cache and
the public upstream density/velocity/Parker field for an exactly equivalent
equatorial geometry.  `1D3D02` compares the complete MHD shock state in the
spherical +X limit, including shock existence, Mach number, compression,
upstream/downstream vectors, density, and pressure.  Both must agree to
roundoff-level tolerances.

The common-core refactor has since been extended by the shared region,
acceleration, and SEP-interface layers.  `1D3D03` now validates the shared
SOURCE acceleration record, while `SEP03` validates the complete AMPS-facing
1-D/3-D `SEPSourceState` serialization in the same spherical/+X limit.

## Validation

Build the validation executable from `test/`:

```sh
cd test
make
./output/test_swcme --list
./output/test_swcme --all
```

The finite shock geometry is covered by `GEO01`-`GEO08` in
`test/3d/test_geometry.cpp`.  See `test/README.md` for detailed test purposes,
reference calculations, and acceptance criteria.

## Ideal-MHD shock existence and downstream state

The production model now separates the existence of a geometric CME front from
the existence of a physical fast shock.  A local surface point is classified as
a shock only when its normal speed exceeds the upstream normal flow by more
than the local oblique fast-mode speed.  If that criterion is not met, the
model returns `has_shock=false`, `compression=1`, and an unchanged downstream
state.  The legacy `sheath_comp_floor` parameter no longer changes the physical
shock compression and therefore cannot manufacture a shock.

The shared `swcme_shock.hpp` solver evaluates the ideal-MHD Rankine-Hugoniot
conditions in the shock frame.  For a trial density compression it enforces
mass conservation, tangential momentum conservation, and tangential electric
field continuity; normal momentum gives the downstream pressure and a
bracketed scalar solve enforces total-energy-flux conservation.  Accepted
solutions must be compressive, have positive downstream pressure, increase the
entropy proxy, and satisfy the stored conservation residual tolerances.

`swcme3d::Model::shock_state_direction()` is the preferred 3-D API for local
shock diagnostics.  It returns the shock-surface radius and normal, normal shock
speed, fast Mach number, `theta_Bn`, density compression, complete upstream and
downstream primitive states, and conservation residuals.  Upstream quantities
are always evaluated at the actual shock surface rather than at an arbitrary
query radius.

### Surface-owned shock state and query-point invariance

The shock state is now explicitly **owned by the physical shock surface**.  The
Cartesian field evaluators, connectivity solver, directional diagnostics, and
shock-mesh builder all obtain local shock properties through the same
`shock_state_direction()` path.  In particular, ambient density and magnetic
field used to form the Mach number are sampled at `R_shock * u`, not at the
radius of a point where a caller happens to request `n`, `V`, or `B`.  This
prevents the same physical shock from acquiring different compression ratios
when queried from different upstream/downstream locations.

`local_oblique_rc()` remains only as a source-compatible scalar wrapper for
older callers.  Its historical `r_eval_m`, `Rdir_m`, and `n_hat` inputs no
longer control shock physics; the wrapper recomputes the canonical surface state
from the direction and returns scalar projections of that state.  New production
code should call `shock_state_direction()` directly.

The mesh builder and `diagnose_direction()` were also converted to consume the
canonical `LocalShockState` directly, eliminating the last internal secondary
shock-strength path.  Validation tests `SHK13` and `SHK14` permanently guard
query-radius invariance and mesh/diagnostic consistency.

The canonical 3-D shock API retains the exact MHD upstream/downstream state at
the mathematical shock surface.  Transport-facing FULL_ICME fields may represent
that discontinuity with the finite C1 shock layer selected by
`RESOLVED_COMPRESSION`; the inner edge of that numerical layer is pinned to the
exact RH downstream state.  SOURCE mode instead uses SHOCK_ONLY and therefore
contains no resolved RH compression in the transport background.  The interior
sheath relaxation remains phenomenological and separate from the exact shock
diagnostic state.

The 1-D model uses the same shared ideal-MHD jump solver.  Its radial direction
is the shock normal and the Parker azimuthal field is tangential to that normal.
This removes the former 1-D compression-floor shock and makes the downstream
normal speed satisfy the same physical jump conditions as 3-D.

The shock validation suite `SHK01`-`SHK14` covers shock/no-shock classification,
obliquity, independent parallel and perpendicular limiting solutions, oblique
branch continuity, mass flux, normal magnetic field, tangential electric field,
momentum and energy fluxes, entropy/admissibility, the weak-shock limit, query-radius invariance, and canonical-state consistency across diagnostics/mesh outputs.

## Correct shock-surface mesh topology and area-weighted sampling

`Model::build_shock_mesh()` now constructs an explicit triangular manifold
instead of triangulating a rectangular theta-phi array.  The old layout stored
`phi=0` and `phi=2*pi` as separate vertices and stored an entire azimuthal ring
at `theta=0`, even though every point in that ring is the same physical apex.
Those duplicate nodes necessarily generated zero-area apex triangles and a
duplicated periodic seam.  Removing zero-area cells after construction would
not repair the topology because adjacency and surface integration would still
contain duplicate physical vertices.

The corrected topology is shape aware:

- a finite SSE cap contains one unique apex plus `nTheta` rings of exactly
  `nPhi` unique vertices; the final ring is the physical half-width boundary;
- a Sphere or Ellipsoid contains one unique apex/north pole, `nTheta-1`
  periodic interior rings, and one unique rear/south pole;
- no ring contains a duplicate `phi=2*pi` point.  Periodicity is represented
  only by triangle connectivity using `(iphi+1)%nPhi`;
- every triangle is wound outward.  `compute_triangle_metrics()` rejects
  repeated-index, non-finite, numerically degenerate, or inward-oriented cells
  rather than returning zero area or a fabricated normal.

`ShockMesh` records `n_theta_intervals`, `n_phi`, and `closed_surface` so a
consumer can audit the topology directly.  Existing 1-based triangle indices
remain unchanged for Tecplot compatibility.

Area weighting is now a production API rather than a demo-side convention.
`build_area_sampling_table()` forms a strictly monotone cumulative distribution
from physical triangle areas using long-double accumulation, and
`sample_triangle_by_area(table,u)` maps a caller-supplied `u` in `[0,1)` to a
triangle.  SWCME intentionally does not own the random-number generator: AMPS
or another caller controls the RNG/seed, while SWCME guarantees that a spatially
uniform source is sampled in proportion to physical surface area rather than
angular-grid index.  Zero/non-finite cell areas are rejected and never skipped.

The deterministic validation gates `MSH01`-`MSH05` cover nondegeneracy, outward
orientation, second-order surface-area convergence, unique apex/periodic-seam
topology, and fixed-seed chi-square tests of area-weighted stochastic patch
selection.  See `MESH_TOPOLOGY_FIX_NOTES.md` for construction formulas, node/
cell counts, and the validation rationale.

## Correct velocity-divergence treatment

Velocity divergence is now evaluated with an operator that matches the actual
velocity field rather than applying one radial formula to every 3-D state.  The
production implementation is centralized in `swcme_divergence.hpp` and the
analytical radial-profile derivatives live beside the region interpolation in
`swcme_regions.hpp`.

For the controlled `SHOCK_ONLY` background,

```text
V = V_sw e_r,  |V_sw| = constant,
div(V) = 2 V_sw / r.
```

Both 1-D and 3-D use this exact expression.  No finite-difference stencil is
used, so the SEP adiabatic energy-change term has no step-size noise in the
baseline Parker/Leblanc wind.

The 1-D `FULL_ICME` velocity remains purely radial.  Its divergence is therefore
computed analytically as

```text
div(V) = 2 V_r/r + dV_r/dr.
```

`dV_r/dr` is the exact derivative of the same shock/sheath/LE/TE smoothstep
profile that returns `V_r`; value and derivative are generated by one common
`RadialVelocityState`.  The legacy `dr_frac` argument remains in the 1-D public
signature for source compatibility but is not used.

The 3-D `FULL_ICME` field is different.  An oblique Rankine-Hugoniot jump can
produce tangential velocity components, and finite SSE/ellipsoid surfaces make
the region state vary with direction.  A ray-only expression
`r^-2 d(r^2 V_r)/dr` therefore omits physical angular/tangential contributions.
The canonical 3-D path now evaluates

```text
div(V) = dVx/dx + dVy/dy + dVz/dz
```

with a shared second-order Cartesian operator.  Centered differences are used
when the stencil lies in the model domain; if one side crosses the explicit
`1.05 R_sun` inner boundary, a documented second-order one-sided stencil is
used instead of clipping the radius.

The preferred 3-D APIs are:

```cpp
swcme3d::Model::compute_divV_checked(...)
swcme3d::Model::compute_divV_cartesian_checked(...)
```

`compute_divV_checked()` selects the exact analytical result in `SHOCK_ONLY`
and the full Cartesian result in `FULL_ICME`.
`compute_divV_cartesian_checked()` always evaluates the numerical full-vector
operator and is useful for convergence studies.  The old
`compute_divV_radial_checked()` symbol is retained for source compatibility but
now delegates to the corrected mode-aware implementation; it no longer applies
a radial approximation to non-radial `FULL_ICME` flow.

The validation gates are `DIV01`-`DIV03`: exact constant-radial divergence, a
manufactured nonconstant radial profile with `<1e-10` relative error, and
second-order convergence of the general Cartesian operator.  See
`VELOCITY_DIVERGENCE_FIX_NOTES.md` and `test/README.md` for the fixtures and
acceptance criteria.

## Completed remediation baseline

The planned SWCME remediation baseline now includes common constants/units and
configuration validation, Parker/Leblanc ambient physics, shared kinematics,
finite shock geometry, the ideal-MHD jump, cobpoint connectivity, common
sheath/ejecta handling, a mutually exclusive acceleration representation,
explicit numerical status propagation, corrected shock-mesh topology,
validated velocity divergence, standardized defaults/scope, and the
SWCME-to-SEP integration/campaign layer described above.  Further development
can therefore be treated as new physics or mission/application integration
rather than completion of the original numerical-remediation list.

## Repaired sheath/ejecta region model and SHOCK_ONLY/FULL_ICME modes

The phenomenological downstream-region model is now shared through
`swcme_regions.hpp`.  This closes the former 1-D/3-D divergence in layer
geometry, ejecta factors, and artificial leading/trailing-edge smoothing.

Two explicit modes are available through `Params::region_mode`:

- `swcme::regions::Mode::ShockOnly` returns the undisturbed analytical
  Parker/Leblanc solar-wind state everywhere.  Shock geometry, connectivity,
  and local MHD shock/source diagnostics remain available through their
  dedicated APIs, but no sheath/ejecta plasma modification is applied to the
  transport-facing background.  This is the recommended controlled baseline
  for the SEP connectivity/perpendicular-diffusion study.
- `swcme::regions::Mode::FullICME` adds the optional phenomenological sheath and
  magnetic-ejecta profile behind the geometric CME/shock surface.

The exact physical shock remains an explicit surface in the diagnostic API and
its downstream state is the ideal-MHD Rankine-Hugoniot solution.  For transport,
FULL_ICME is paired with `RESOLVED_COMPRESSION`: a symmetric finite-width C1
layer is centered on the mathematical shock, reaches the analytical upstream
state at its outer edge, and reaches the exact RH downstream state at its inner
edge.  The sheath then relaxes smoothly toward a leading-edge target.  The
empirical `sheath_comp_floor` no longer participates in shock or region physics
and is retained only for source compatibility.

### Self-similar local layer geometry

The public sheath/ejecta thickness inputs are specified as AU at a 1-AU shock.
They are now interpreted as dimensionless self-similar fractions.  For a local
shock-surface radius `R_sh(u)`,

```text
R_LE = (1 - f_sheath) R_sh
R_TE = (1 - f_sheath - f_ejecta) R_sh.
```

This is applied to the **local** SSE/ellipsoid radius, not the apex radius.
Consequently the shock, leading-edge, and trailing-edge surfaces remain nested
and geometrically similar from apex to flank.  Configurations with
`f_sheath + f_ejecta >= 1` are rejected instead of repaired by runtime clipping.

The LE and TE smoothing inputs are likewise local self-similar fractions. Each
configured width is the total width of a symmetric C1 smoothstep transition
centered on the nominal boundary. Oversized widths are rejected before
preparation; accepted widths are never capped or otherwise changed.

### Sheath and magnetic-ejecta targets

For a physical fast shock, the sheath starts at the complete MHD downstream
state and relaxes toward the local Parker/upstream state at the leading edge.
If a geometric CME front exists but the fast shock has decayed, the model does
not fabricate a sheath compression; that portion of the profile remains
ambient while the optional ejecta may still be represented.

Magnetic-ejecta factors are now honored exactly:

```text
n_ME = f_ME * n_up
V_ME = V_ME_factor * V_sw.
```

Values below unity therefore correctly produce density depletion and slower
bulk ejecta.  Negative factors are rejected by centralized configuration
validation rather than clipped in the evaluator.  The baseline ejecta magnetic
field remains Parker; a flux-rope/ejecta-field model is intentionally outside
the controlled one-year scope.

Artificial LE and TE interfaces are C1 smooth.  The numerical shock layer is
also C1, but only when `RESOLVED_COMPRESSION` is selected.  Its total width is
`edge_smooth_shock_AU_at1AU * R_sh(local)`; configuration validation limits the
request so it cannot consume the entire sheath. In SOURCE mode that shock-layer
width is forced to zero and the SHOCK_ONLY background remains analytical on
both sides of the mathematical source surface.

Validation tests `REG01`-`REG05` cover SHOCK_ONLY identity, the exact RH state at
the inner edge of the resolved shock layer, sub-unity ejecta factors, local
self-similar surface nesting, and continuity/smoothness of the modeled region
transitions.

See `REGION_MODEL_FIX_NOTES.md` and `test/README.md` for detailed equations,
configuration conventions, and validation procedures.


## Single shock-acceleration representation

`swcme_acceleration.hpp` defines one authoritative acceleration switch shared by
1-D and 3-D:

```cpp
swcme::acceleration::Mode::Source
swcme::acceleration::Mode::ResolvedCompression
```

The model deliberately does **not** expose independent `use_dsa_source` and
`use_compression_acceleration` booleans.  A single enum prevents a run from
enabling both representations of first-order shock acceleration for the same
particle population.

### SOURCE mode

`SOURCE` is the recommended baseline for the connectivity/perpendicular-
diffusion experiment.  Centralized validation currently requires

```text
acceleration = SOURCE
regions      = SHOCK_ONLY
```

so the transport-facing velocity remains the analytical solar wind through the
mathematical shock surface.  A physical fast shock instead produces a
`ShockAccelerationState` with `source_enabled=true`, local shock position and
normal, normal shock speed, compression, `theta_Bn`, fast Mach number, upstream
density and magnetic-field magnitude, and the diagnostic DSA phase-space slope

```text
q = 3 r_c / (r_c - 1).
```

The baseline `relative_source_weight_per_area` remains dimensionless and is
intended only for controlled relative weighting.  The AMPS-facing
`swcme_sep_source.hpp` layer now makes normalization explicit: its default
`RELATIVE_ONLY` mode preserves this dimensionless convention, while optional
`REFERENCE_DIFFERENTIAL_INTENSITY` requires a physical `J(E_ref)` with declared
SI units.

### RESOLVED_COMPRESSION mode

Centralized validation currently requires

```text
acceleration = RESOLVED_COMPRESSION
regions      = FULL_ICME
edge_smooth_shock_AU_at1AU > 0
```

and no prescribed DSA source is exposed.  The exact RH discontinuity remains in
`shock_state_direction()` / the 1-D cached jump for diagnostics, while the
transport-facing primitive fields use one common symmetric C1 shock profile.
For a total local width `w_sh`, the profile spans

```text
R_sh + w_sh/2   upstream endpoint
R_sh            midpoint
R_sh - w_sh/2   exact RH downstream endpoint.
```

The smoothstep derivative vanishes at both endpoints.  The phenomenological
sheath begins at the inner endpoint, and its own profile also has zero slope
there, so there is no second numerical kink at the resolved-shock/sheath join.
The same `swcme_regions.hpp` formulas are called by 1-D and 3-D, including at
finite-SSE/ellipsoidal flanks where `w_sh` scales with the **local** shock radius.

In `RESOLVED_COMPRESSION`, `ShockAccelerationState::source_enabled` is false,
`resolved_compression_enabled` is true, the active DSA slope is intentionally
reported as unavailable, and the source weight is zero.  This prevents a
consumer from silently applying a pre-imposed DSA source on top of the resolved
`div(V)` accelerator.

The deterministic helper `swcme::acceleration::serialize_csv()` exists for
regression/audit comparisons.  It is not the final science-data format.  Tests
`ACC01`-`ACC05` and `1D3D03` verify the no-double-counting gate, C1 resolved shock
profile, RH endpoints, 1-D/3-D smoothing identity, disabled source in resolved
mode, and byte-identical SOURCE records in the spherical radial limit.

See `SHOCK_ACCELERATION_FIX_NOTES.md` and `test/README.md` for the validation
fixtures and implementation details.

## Explicit numerical status and error propagation

`swcme_status.hpp` provides the shared `swcme::ModelStatus` contract used by
numerically sensitive 1-D/3-D evaluation paths.  The production physics layer
no longer converts an invalid intermediate state into zero, ambient plasma, a
clipped radius, or an arbitrary +X direction merely to keep a calculation
running.

New code should prefer the checked APIs:

```cpp
swcme1d::Model::evaluate_radii_fast_checked(...)
swcme1d::Model::evaluate_radii_with_B_div_checked(...)

swcme3d::Model::shock_state_direction_checked(...)
swcme3d::Model::evaluate_cartesian_fast_checked(...)
swcme3d::Model::evaluate_cartesian_with_B_checked(...)
swcme3d::Model::evaluate_cartesian_with_B_div_checked(...)
swcme3d::Model::compute_divV_checked(...)
swcme3d::Model::compute_divV_cartesian_checked(...)
// compute_divV_radial_checked(...) remains as a compatibility alias
```

They return `ModelStatus` with an explicit code, context, and failed sample
index for batch queries.  The existing void/bool science wrappers remain for
source compatibility but convert numerical failures to exceptions instead of
silently repairing the result.  `NO_SURFACE` remains an expected geometrical
outcome for directions outside a finite front and is distinguishable from an
actual numerical failure.

The analytical Leblanc/Parker lower boundary, `1.05 R_sun`, is now a real model
domain: science queries below it return `OUTSIDE_MODEL_DOMAIN`.  The low-level
solar-wind equations no longer clip the caller's radius.  Vector norms use
`std::hypot` and degenerate/non-finite vectors are rejected; the previous +X
normalization fallback has been removed.

The shared ideal-MHD jump result also carries an explicit `SolveStatus`
(`NO_SHOCK`, `SOLVED`, `INVALID_INPUT`, `NO_PHYSICAL_BRACKET`,
`INVALID_ACCEPTED_STATE`, or `CONSERVATION_FAILURE`).  A failed super-fast RH
solve is propagated as `SHOCK_SOLVER_FAILURE`; it is not converted to a
compression-one ambient state.

Mesh and checked Tecplot output paths reject non-finite physics instead of
writing a sanitized surrogate.  `ERR01`-`ERR05` protect outside-domain behavior,
non-finite batch input, explicit RH outcomes, degenerate-vector rejection, and
writer rejection of corrupt mesh data.  See `NUMERICAL_STATUS_FIX_NOTES.md` and
`test/README.md` for the complete status semantics and validation fixtures.
