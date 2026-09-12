# SWCME validation

The validation suite is a standalone C++ executable that calls production
SWCME interfaces. It uses the repository's existing Make-based build approach;
no second build system or external test dependency is required.

## Directory layout

```text
test/
  core/       Common registry, runner, reporting, CFG/DEN/KIN/SHK tests
  1d/         Tests specific to the production 1-D model
  3d/         Tests specific to the production 3-D model
  reference/  Reviewed reference/event data used by comparison tests
  profiles/   SMOKE/ROUTINE/FULL/EVENT test selections
  python/     Regression tests for the campaign manager
  output/     Generated executables and campaign artifacts (ignored by Git)
  run_tests.py  Python campaign manager
  event_config.example.json  executable EVENT/sweep example
```

Every validation is linked into the single `output/test_swcme` executable.

## Build

From `srcSEP/swcme`:

```sh
make -C test clean all
```

From `srcSEP/swcme/test`, the equivalent command is `make clean all`.

To build the test executable and demonstrations in parallel and then run the
complete registered validation suite, run the following command from
`srcSEP/swcme/test`:

```sh
make -j test
```

The `test` target first brings all required binaries up to date and then runs
`output/test_swcme` in registry order. Because `-j` does not set an explicit
job limit, use `make -jN test` instead when the build host should be limited to
`N` concurrent compilation jobs.

## Command-line interface

```sh
./output/test_swcme                 # run all tests; same as --all
./output/test_swcme --all           # run all tests in registry order
./output/test_swcme --list          # list tests without executing them
./output/test_swcme --test PST01    # prepared-state immutability
./output/test_swcme --test PST02    # prepared-state ownership rejection
./output/test_swcme --test PST03    # configuration-state ownership rejection
./output/test_swcme --test PST06    # prepared-record integrity rejection
./output/test_swcme --test PST04    # concurrent prepared-state evaluation
./output/test_swcme --test OUT02    # checked output-failure propagation
./output/test_swcme --test OUT03    # transactional output commit
./output/test_swcme --test OUT05    # model-domain output preflight
./output/test_swcme --test OUT04    # BoxSpec structural validation
./output/test_swcme --test OUT06    # mesh and metric output validation
./output/test_swcme --test OUT01    # independent Tecplot parsing
./output/test_swcme --test OUT07    # build/run/parse all demonstrations
./output/test_swcme --test PST07    # direct/AMPS adapter equivalence
./output/test_swcme --test PST05    # prepared-state lifetime and relocation
./output/test_swcme --test PST08    # state-ownership performance guardrails
./output/test_swcme --test OUT08    # strict-writer record semantics
./output/test_swcme --test CFG03    # smoothing-width rejection policy
./output/test_swcme --test CFG04    # configured-radius domain policy
./output/test_swcme --test KIN09    # kinematic extrapolation domain
./output/test_swcme --test CFG01    # run exactly CFG01
./output/test_swcme --test CFG02    # run exactly CFG02
./output/test_swcme --test DEN01    # run exactly DEN01
./output/test_swcme --help          # usage, options, and examples
```

Options are mutually exclusive. Unknown test IDs, unknown options, a missing
value after `--test`, and combined modes return exit code 2. A completed test
run returns 0 when every mandatory test passes and 1 when any test fails.
Subcheck SKIPs identify unavailable production paths and do not by themselves
fail an otherwise successful test.

The registry in `core/test_swcme.cpp` is the only source for `--list`, `--all`,
and test lookup, so displayed and executed tests cannot silently diverge.

## PST01: prepared-state immutability

`PST01` verifies that a successfully prepared 1-D or 3-D state is independent
of every later attempt to modify its model configuration.  The model remains
configurable during setup, but its first successful `prepare_step()` freezes
the configuration.  Failed preparation leaves the setup phase unlocked so an
invalid model can be corrected and retried.

The 1-D fixture evaluates density, velocity, Parker-field components, field
magnitude, and divergence at three fixed radii and serializes all 18 values.
It then attempts every retained legacy mutation path:

- `SetParams`, `SetCME`, and both kinematics setters;
- `SetAmbient`, including the historically defective `sin_theta` case;
- region and shock-acceleration mode setters;
- geometry, smoothing, and sheath/ejecta setters;
- copy assignment.

Each attempt must throw `std::logic_error` with the immutable-lifecycle
diagnostic before changing the configuration digest or model identity.  The
same state is reevaluated after every attempt, and its serialized result must
be bitwise identical to the baseline.

Compile-time API inspection separately requires the legacy 1-D
`MutableParams()` raw-reference escape hatch to be unavailable.  Merely checking
that accessor at call time would be insufficient because a reference obtained
before preparation could be retained and used after the model freezes.

The 3-D fixture uses compile-time inspection to confirm that no mutable Params
accessor exists, then verifies that copy assignment is rejected after
preparation.  Density, all velocity and magnetic-field components, and
divergence at three Cartesian points remain bitwise identical.  Both fixtures
exercise `reconfigured(params)`: it must create a new owner and configuration
snapshot, permit independent preparation and different physics, and leave the
original state unchanged.

Run the gate directly with:

```sh
./output/test_swcme --test PST01
```

`PST01` is the first test in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST02: cross-model prepared-state rejection

`PST02` verifies that each prepared state belongs to the exact model instance
that created it.  This is intentionally stronger than configuration equality:
two separately constructed models with identical parameters receive different
process-local identities, and neither accepts the other's `StepState`.

The fixture prepares a state with Model A and attempts to consume it with Model
B using equal configurations and with Model C using a different solar-wind
configuration.  It covers:

- checked 1-D background, full-field, shock-source, and Tecplot-writer paths;
- checked 3-D background, magnetic-field, divergence, shock, acceleration, and
  Tecplot bundle paths;
- legacy 3-D geometry, mesh, and connectivity wrappers; and
- 1-D/3-D AMPS-facing background, point source, surface source, and cobpoint
  source adapters.

Every checked call must return `STATE_MODEL_MISMATCH` with
`has_model_identities=true`, the receiver's identity in
`expected_model_identity`, and Model A's identity in
`supplied_model_identity`.  Numerical arrays and source/background objects are
initialized with sentinels and must remain unchanged.  Pre-existing 1-D and
3-D output files contain fixed byte strings and must remain byte-for-byte
identical, proving rejection occurs before `fopen("w")`.  Legacy wrappers must
throw an exception whose diagnostic retains `STATE_MODEL_MISMATCH`; they may
not turn ownership misuse into a valid radius, mesh, or disconnected result.

Run the gate directly with:

```sh
./output/test_swcme --test PST02
```

`PST02` is included in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.  A default-constructed state has owner identity zero
and is rejected by the same contract.

## PST03: cross-configuration state rejection

`PST03` verifies that prepared-state provenance protects the complete resolved
configuration, not only the `gamma_ad` value that originally exposed the
defect.  `prepare_step()` records an allocation-free deterministic digest in
the state.  Every state-consuming API compares the receiver's current digest
before it evaluates physics or modifies caller-owned output.

The test prepares one reviewed 3-D baseline, then constructs receivers that
differ by exactly one field representing each validation-plan family:

- adiabatic index (`gamma_ad`);
- geometry (`axis_ratio_y`, including inactive-field coverage);
- Parker orientation (`solar_rotation_axis`);
- kinematics mode;
- thermal-closure input (`T_K`);
- Parker field/polarity normalization input (`B1AU_nT`); and
- region mode.

Parker radial polarity and the proton-only pressure closure are compile-time
resolved conventions rather than independently mutable `Params` fields.  They
are therefore included as explicit versioned digest tags; `B1AU_nT` and `T_K`
exercise their runtime normalization/closure inputs.  A multi-field receiver
checks that no parameter-specific rejection branch exists.

For every foreign receiver, the checked evaluator must return
`STATE_MODEL_MISMATCH` before writing its sentinel outputs.  PST02 ownership has
intentional precedence, while `has_configuration_digests=true` and the
receiver/prepared digest pair prove the additional configuration mismatch.
An explicitly populated default-equivalent receiver must have the same digest
as the baseline yet remain rejected as a foreign model.

Because PST01 rejects supported post-prepare model mutation, PST03 deliberately
alters the digest tag only on a copied negative-test state.  The same-owner
consumer must return `STATE_CONFIGURATION_MISMATCH`, carry current and supplied
digests, and leave density/velocity sentinels unchanged.  A reviewed golden
digest plus an independently constructed equal configuration guards
reproducibility across runs and compiler rebuilds.

Run the gate directly with:

```sh
./output/test_swcme --test PST03
```

`PST03` follows `PST02` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST06: prepared-state record integrity

`PST06` verifies that model ownership and configuration equality cannot be
bypassed by modifying a cached `StepState` value.  Each successful preparation
stores a private seal produced by explicit field serialization.  The seal
covers canonical common solar-wind/kinematic state, region boundaries,
acceleration configuration, 1-D shock/RH data, 3-D frame and geometry caches,
and every top-level legacy compatibility mirror.

Compile-time assertions require `integrity_digest()` to return by value rather
than expose a mutable reference.  They also require both state types to remain
copy- and move-constructible.  Runtime fixtures exercise copy construction,
move construction, copy assignment, and move assignment; valid copies retain
the original owner, seal, and evaluability.

The negative matrix copies a valid state and changes exactly one record field
at a time.  It includes every top-level 1-D and 3-D mirror plus representative
members of every nested canonical record.  Each checked evaluator must return
`STALE_PREPARED_STATE`, set `has_state_integrity`, report the prepared and
recomputed seals, and leave all numerical output sentinels unchanged.  A legacy
1-D wrapper must throw with the same status name.  Separate owner/configuration
tag tests remain in PST02/PST03 because those higher-precedence diagnostics are
intentional.

Run the gate directly with:

```sh
./output/test_swcme --test PST06
```

`PST06` follows PST03 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST04: concurrent prepared-state evaluation

`PST04` verifies the production usage pattern required by AMPS: configuration
and `prepare_step()` finish on one setup thread, after which multiple workers
read the same model and the same immutable prepared state.  It does not grant
permission to call model mutators, prepare a new state concurrently, share
caller-owned destination arrays, or invoke file writers from multiple threads.

One 1-D and one 3-D state are each prepared once.  The fixture then exercises
nine operation families:

- scalar 1-D and 3-D AMPS background queries;
- direct 1-D and 3-D full-field batch evaluators;
- 1-D shock-source conversion;
- 3-D directional shock and source evaluation;
- 3-D Parker-line connectivity with complete cobpoint shock records; and
- simultaneous, intentionally different 1-D/3-D failures that retain their
  own status code, context, sample index, offending-value flag, and outputs.

Every public result is serialized field by field into fixed-width words.  The
serializer includes every `ModelStatus` field, all background/source/shock
numbers and flags, every connectivity root, and all batch output elements.  It
does not compare raw structure memory, so unspecified padding cannot cause a
false race diagnosis.  Each concurrent result must be bitwise identical to a
serial oracle; no roundoff allowance is currently required.

The workload runs eight repetitions at 1, 2, 4, and 8 threads under
forward-interleaved, reverse-interleaved, and operation-grouped job orders.  A
one-shot start gate creates overlap, while each output slot has exactly one
writer.  The test framework is used only after joining the workers, preventing
the harness itself from introducing a race.  Finally, both prepared-state
integrity seals are recomputed to prove the stress run did not mutate them.

Run the normal concurrency gate with:

```sh
./output/test_swcme --test PST04
```

When the compiler and runtime support ThreadSanitizer, rebuild the entire
executable with instrumentation so both the test and production implementation
are observed:

```sh
make clean
make CXXFLAGS="-O1 -g -std=c++17 -Wall -Wextra -Wpedantic -fsanitize=thread -fno-omit-frame-pointer"
TSAN_OPTIONS=halt_on_error=1 ./output/test_swcme --test PST04
```

A supported ThreadSanitizer run must report no race.  Some container kernels
cannot initialize ThreadSanitizer's shadow memory; that platform limitation is
not a passing sanitizer result and should be recorded as unavailable.  `PST04`
follows `PST06` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## OUT02: write failure detection and propagation

### What is tested

`OUT02` verifies the complete lifecycle of every production Tecplot product:
the 1-D radial profile and shock-versus-time history, plus the 3-D surface,
four-zone surface/volume/face bundle, and standalone min-X face.  It
distinguishes failure to open a destination from failure after opening, checks
exact partial-write byte offsets, exercises flush, persistent stream-error,
and close failures, and verifies that cleanup does not overwrite the first
diagnostic.  Both the inline/header 1-D implementation and the separately
compiled 3-D implementation are exercised.

### Why it is tested

Buffered output can accept formatted data into memory and report the actual
device failure only at `fflush()` or `fclose()`.  The former writers ignored
formatted-write and close results, so a full device or interrupted write could
produce a truncated science file while returning success.  That makes a
partially written zone indistinguishable from a complete model product to an
automated observational campaign.

### How it is tested

The deterministic portion supplies a C-compatible `FileOperations` table to
the checked writer APIs.  An isolated in-memory `FaultSink` first verifies a
successful open/write/flush/error/close sequence, then injects:

- failure to open;
- rejection of the first formatted write;
- a partial write after exactly 57 bytes, combined with a later close failure;
- a separately compiled 3-D partial write after exactly 41 bytes;
- failure at flush;
- a persistent stream error after a successful flush; and
- failure at close.

Each post-open case must return `FILE_WRITE_FAILURE`, set
`has_io_byte_offset`, report the exact first unaccepted byte, identify the
failed record or lifecycle phase in `context`, and close the acquired handle.
The combined partial-write/close case proves that cleanup cannot mask the
earlier failure.

When `/dev/full` is available, the shared direct-stream layer is exercised
through the default `stdio` backend and must report `FILE_WRITE_FAILURE` with
byte context.  The probe is below OUT03's regular-file transaction boundary:
model writers stage a sibling file and never replace a device.  A missing
`/dev/full` is recorded as one explicit platform skip rather than silently
weakening the deterministic fault matrix.

### Expected result

The test passes only when open rejection is `FILE_OPEN_FAILURE`, every
post-open failure is `FILE_WRITE_FAILURE` and never `OK`, byte/row/zone context
identifies the first failure, and all acquired handles are closed.  Normal
injected 1-D and 3-D writes must remain successful.  OUT02 detects partial
output; OUT03 independently proves that such staged output is never published.

Run the gate directly with:

```sh
./output/test_swcme --test OUT02
```

`OUT02` follows `PST04` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT03: transactional output commit

### What is tested

`OUT03` verifies that all model products are written under a private,
same-directory staging name and become visible at the requested path only
after write, flush, stream-error, and close checks succeed.  It covers
exclusive staging creation, collision retry, successful replacement of an
existing regular file, creation of a new file, refusal to replace a
non-regular object, failed atomic commit, staging cleanup, and preservation of
the first failure when cleanup itself fails.

### Why it is tested

OUT02 prevents a truncated file from being labeled successful, but detection
alone is insufficient when a writer opens the final path directly: the old
complete product has already been truncated.  An observational campaign may
then find a corrupt file even though the status correctly reports failure.
OUT03 makes publication transactional so downstream readers see the previous
complete result until a complete replacement is ready.

### How it is tested

The deterministic `FaultSink` models separate staging and destination byte
stores.  Successful injected transactions are run for the 1-D radial profile,
1-D shock history, 3-D surface, 3-D bundle, and 3-D face.  The matrix also
verifies that:

- the staging path begins with the exact destination plus `.swcme-tmp-` and is
  opened with exclusive mode `wx`;
- one simulated staging-name collision causes a retry and then succeeds;
- successful close is followed by commit, with no failure cleanup;
- partial-write and close failures skip commit, remove staging, and preserve
  the destination sentinel byte-for-byte;
- injected commit failure returns `FILE_COMMIT_FAILURE`, records the complete
  staged byte count, preserves the destination, and removes staging;
- staging-open failure never invokes commit or removal; and
- failed cleanup cannot mask the earlier `FILE_WRITE_FAILURE`.

Production-filesystem integration then replaces an existing sentinel file,
installs a brand-new 3-D output file, and verifies that neither success leaves
a `.swcme-tmp-*` sibling.  A legacy boolean writer must perform the same
complete replacement.  A directory is used as a safe non-regular destination:
the checked API must return `FILE_COMMIT_FAILURE`, the legacy API must return
`false`, both must preserve the directory, and neither may leave staging.

### Expected result

Every successful checked or legacy writer exposes one complete final file and
no staging file.  Any pre-commit failure leaves the previous destination
unchanged and attempts staging cleanup.  Non-regular targets are never
replaced.  Status codes distinguish staging-open, stream, and commit failures,
and later cleanup cannot overwrite the first diagnostic.  The guarantee is
atomic visibility on supported POSIX filesystems; crash durability through
file and directory `fsync()` is outside OUT03.

Run the gate directly with:

```sh
./output/test_swcme --test OUT03
```

`OUT03` follows `OUT02` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT05: model-domain output preflight

### What is tested

`OUT05` verifies that every public Tecplot writer validates its complete
requested model domain before any output backend callback or filesystem action.
Coverage includes a non-first invalid 1-D radius, a non-finite 1-D coordinate,
a non-finite precomputed 1-D field, a late out-of-range data-driven history
time, a non-first invalid shock-surface vertex, an interior invalid structured
volume point, and an invalid min-X face.  It also checks the inclusive
lower-radius boundary and the legacy boolean face wrapper.  OUT04 separately
owns malformed `BoxSpec` structure and arithmetic-range validation.

### Why it is tested

OUT02 reports incomplete writes and OUT03 prevents incomplete staging data from
replacing a destination, but neither guarantee by itself proves that invalid
physics requests are side-effect free.  Preflight is important to automated
campaigns because a late bad sample should be rejected before creating files,
invoking a custom backend, or performing expensive partial output evaluation.
The originating model-domain status and row must remain available so the bad
campaign sample can be corrected without interpreting an I/O error.

### How it is tested

The deterministic `FaultSink` counts open, commit, and remove callbacks while
holding a sentinel destination.  Each invalid request is passed to the same
checked API used by production and must return with `open_attempts == 0`, no
commit, and no cleanup callback.  Assertions verify the status class,
`sample_index`, offending value, and zone-specific context:

- radial coordinates distinguish `NONFINITE_INPUT` from
  `OUTSIDE_MODEL_DOMAIN`, while supplied field data use `NONFINITE_RESULT`;
- the data-driven history accepts its first two times but reports the third,
  out-of-range request before opening;
- the surface reports the altered mesh-node index;
- the bundle reports the flattened interior volume row, demonstrating that the
  entire K/J/I grid—not only endpoints—is scanned;
- the face reports its first invalid effective-grid row, and its legacy wrapper
  returns `false` under the same preflight;
- a face exactly at `MIN_RADIUS_M` reaches open and commits successfully.

A production-backend case starts with an existing regular sentinel file,
submits the invalid middle-radius profile, and then verifies that the sentinel
is byte-for-byte unchanged and no `.swcme-tmp-*` sibling exists.

### Expected result

All invalid datasets fail before the first output open.  Each checked call
returns the domain/result status and first failing row described above; legacy
output returns `false`; the prior destination remains unchanged; and no staging
file is created.  The exact lower domain boundary remains writable and commits
one complete product, proving the preflight does not reject valid edge points.

Run the gate directly with:

```sh
./output/test_swcme --test OUT05
```

`OUT05` follows `OUT03` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT04: box specification validation

### What is tested

`OUT04` verifies the shared structural contract for every `BoxSpec` consumer.
The matrix covers non-finite centers and half extents, a negative half extent,
each grid dimension below two, overflowed center-plus/minus-extent bounds, a
doubled span that overflows despite finite bounds, and an `Ni*Nj*Nk` product
that cannot fit `size_t`.  It exercises both the four-zone bundle and
standalone min-X face, the legacy boolean face wrapper, and the
`default_apex_box()` factory.  A collapsed box with zero half extents is the
positive boundary case.

### Why it is tested

Before OUT04, the documented `Ni,Nj,Nk >= 2` rule was not enforced: the bundle
accepted one sample per axis, and the face did not validate `Ni` at all.
Negative half extents silently reversed grid orientation, while huge finite
values could overflow bounds, the `2*h` coordinate span, or flattened point
cardinality.  Such malformed inputs should be rejected as configuration errors
before OUT05 traverses points or the output layer creates a staging file.

### How it is tested

Invalid boxes are submitted through production checked writer APIs backed by
the deterministic `FaultSink`.  Every case must produce zero open attempts, no
commit or cleanup callback, and leave the simulated destination sentinel
unchanged.  The expected diagnostics are:

- `NONFINITE_INPUT` with an offending value for non-finite members and derived
  bound/span overflow;
- `INVALID_CONFIGURATION` with the offending value for negative extents and
  dimensions below two; and
- `INVALID_CONFIGURATION` for checked point-count multiplication overflow.

The face test sets only `Ni=1` to prove validation is consumer-independent.
The point-count case uses `INT_MAX` on all axes and must return immediately
rather than entering grid traversal.  The legacy writer must return `false`
without creating its requested path.  Factory calls with negative or NaN
half-size and resolution one must throw, while a normal request must return a
three-by-three-by-three box.  A production-backend check verifies an existing
sentinel file remains byte-for-byte unchanged with no `.swcme-tmp-*` sibling.

### Expected result

Every malformed box fails before filesystem or injected-backend access with
the status class described above.  Both checked writers enforce the identical
contract, the legacy wrapper returns `false`, invalid factory inputs throw, and
no destination or staging artifact is created.  Finite nonnegative zero
extents with dimensions of at least two remain valid and commit successfully.

Run the gate directly with:

```sh
./output/test_swcme --test OUT04
```

`OUT04` follows `OUT05` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT06: mesh output validation

### What is tested

`OUT06` verifies that the surface-only and bundle writers accept only a
complete, physically representable, topologically valid, and internally
consistent `ShockMesh`/`TriMetrics` record.  Cases cover empty and mismatched
arrays, non-finite nodal physics, non-unit normals, sub-unity compression,
negative normal speed, zero/out-of-range/repeated connectivity, coincident
vertices, partial metrics, non-finite metrics, finite stale area, and metrics
made stale by a nodal-state change.  It also verifies canonical auto-computation
from completely empty metrics, legacy-wrapper rejection, and production-file
preservation.

### Why it is tested

Previously, writers checked only parallel-array sizes and finiteness.  If
`T.area.size()` already equaled the triangle count, invalid connectivity and
degenerate geometry could bypass `compute_triangle_metrics()` entirely.
Finite but stale cell metrics could likewise be serialized after mesh geometry
or nodal shock values changed.  Such a file is syntactically readable but no
longer represents the supplied surface, which can corrupt area-weighted SEP
source analysis without producing an obvious I/O failure.

### How it is tested

Each invalid record is submitted through the production checked writer with a
deterministic `FaultSink`.  Every rejection must make zero open, commit, and
remove calls and preserve the destination sentinel.  Assertions verify:

- size/empty-record failures return `INVALID_MESH`;
- non-unit normals, `rc < 1`, and `Vsh_n < 0` identify the bad node;
- a non-finite nodal or metric value returns `NONFINITE_RESULT` with its node or
  triangle index;
- a bad non-first connectivity entry returns `INVALID_MESH` with its triangle
  index for zero, out-of-range, and repeated indices;
- coincident coordinates are rejected by canonical triangle-quality checks;
- a partially populated metric record is rejected;
- complete supplied metrics are compared with freshly computed area, normals,
  centroids, mean compression, and mean speed, catching both a modified area
  and an unchanged metric record after nodal `rc` changes;
- an entirely empty metric record triggers canonical computation and commits a
  valid complete surface; and
- the bundle and legacy surface wrapper enforce the same mesh contract.

A production-backend case starts with an existing sentinel file and invalid
connectivity.  The writer must preserve the exact bytes and leave no
`.swcme-tmp-*` sibling.

### Expected result

Every malformed, degenerate, partial, non-finite, or stale mesh record fails
before output access with the documented status and available row index.  No
invalid checked or legacy call creates or changes a destination.  Valid
canonical supplied metrics and completely empty auto-computed metrics both
produce a complete committed surface.

Run the gate directly with:

```sh
./output/test_swcme --test OUT06
```

`OUT06` follows `OUT04` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT01: independent output parsing

### What is tested

`OUT01` validates the complete serialized form of every public 3-D Tecplot
product: the surface-only FETRIANGLE/BLOCK file, the four-zone bundle containing
surface-cell, surface-nodal, structured-volume, and min-X-face data, and the
standalone structured min-X face.  It checks titles, all 25 variable names,
their embedded units and order, zone declarations, variable locations,
structured dimensions, node and element counts, BLOCK lengths, POINT row
widths, numerical finiteness, triangle connectivity, and end-of-file position.

The unit-qualified schema is part of the validated external contract.  Its
coordinate and centroid fields are `[m]`, density is `[m^-3]`, velocity and
shock-speed fields are `[m/s]`, magnetic fields are `[T]`, divergence is
`[s^-1]`, area is `[m^2]`, and compression/normals/reserved directions are
dimensionless `[-]`.

### Why it is tested

A successful writer status proves that bytes were committed, but not that an
independent consumer can interpret those bytes correctly.  A missing zone,
incorrect declared dimension, shifted variable order, omitted unit, wrong row
width, invalid triangle index, nonfinite token, or extra trailing record can
all leave a complete filesystem transaction that is nevertheless unusable or
scientifically misinterpreted.  This test closes that semantic gap before the
demonstration and observational campaigns consume the files.

### How it is tested

The test constructs a deterministic spherical shock mesh and a minimal valid
`2x2x2` apex box, writes all three products through the production stdio and
transactional-commit path, and reads them with a parser implemented solely in
`test/core/tecplot_parser.hpp` and `test/core/test_output_parsing.cpp`.  The
parser uses standard line and numeric conversion primitives and owns its
expected title, variable, unit, zone, and layout declarations.  It does not
reuse `CheckedTextFile`, writer formatting constants, BLOCK emitters, or any
production parsing helper.

For surface zones, every declared node/element count and both copies of the
triangle connectivity are compared with the input mesh.  Representative nodal
coordinates/compression, cell area, and nodal-normal values are compared with
the source records.  For the volume and face, the first Cartesian point is
evaluated directly through the public model API and all eleven meaningful
background fields are compared with the parsed row at `%.9e` precision.  The
standalone face's complete parsed row matrix must equal the face zone embedded
in the bundle.

Four parser-sensitivity mutations are also required to fail: a record appended
after the final zone, a variable unit changed from metres to kilometres, a
zero triangle index, and a 26th value appended to a 25-column POINT row.  These
probes demonstrate that success is not the result of a permissive parser that
ignores precisely the defects OUT01 is meant to find.

### Expected result

All production writes return `OK`.  Each unmodified file parses through its
exact final record with no missing or extra data; the bundle has exactly four
zones; declared sizes match the requested mesh and box; all values are finite;
connectivity is one-based, in range, distinct, and identical to the mesh; and
selected values agree with direct API results within serialization precision.
All four deliberately corrupted byte streams are rejected.

Run the gate directly with:

```sh
./output/test_swcme --test OUT01
```

`OUT01` follows `OUT06` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## OUT07: demonstration program execution

### What is tested

`OUT07` treats `demo1d`, `demo3d_1`, and `demo3d_2` as executable
documentation.  It verifies that all three build with the normal strict-warning
flags, exit successfully from a clean directory, report no stderr diagnostics,
produce exactly the files promised in their comments, and leave no private
transaction file.  It then independently parses every Tecplot and CSV output.

For `demo1d`, the gate checks one 1,200-row, 13-variable radial POINT zone.
Each 3-D example must produce a four-zone bundle with the documented
24-by-48 finite-SSE surface, a 12-by-12-by-12 volume, and a 12-by-12 min-X
face.  The extended example must also produce a six-row predefined point cloud,
ten random samples for each surface triangle, and an ideal-MHD strength table.
Both 3-D examples produce 72-hour plasma and shock histories at five-minute
cadence.

### Why it is tested

Before OUT07, both 3-D examples placed a box face through the solar origin,
outside SWCME's `r >= 1.05 R_sun` domain.  The checked bundle writer therefore
rejected the request, but the examples printed a failure and still returned
zero while claiming the missing file had been written.  The extended example
also spent substantial time constructing a 120-by-240 mesh and ten samples per
large-mesh triangle before reaching that failure.  Its comments described
three output zones although the production bundle has four, and its strength
diagnostic inferred fast Mach number from an obsolete gas-dynamic proxy rather
than using the production ideal-MHD shock state.

Examples are commonly copied into science workflows.  A compile-only test
cannot detect unsupported geometry, ignored output status, missing products,
stale schemas, truncated CSV streams, or comments that no longer match model
physics.

### How it is tested

`make demos` builds the examples as `output/demo1d`, `output/demo3d_1`, and
`output/demo3d_2`; both 3-D programs link the same `output/swcme3d.o` used by
the validation executable.  OUT07 creates a unique run root and one empty
directory per binary, fixes `LC_ALL=C`, captures stdout and stderr, and examines
the process termination code.

The gate compares each directory with an exact expected manifest.  It parses
the radial profile, both bundles, and both point clouds with the OUT01
validation-only parser, whose schema remains independent of production writer
constants.  A separate strict CSV reader verifies exact ordered headers, row
widths, numerical finiteness, final newlines, sample counts, uniform time
spacing, and endpoints.  Surface node/element counts are derived independently
from the documented 24 polar intervals and 48 azimuthal nodes.  Successful
artifacts are removed; failed artifacts and captured logs are retained at
`output/OUT07_demo_runs_<run-id>`.  Process-and-start-time-specific roots allow
simultaneous validation campaigns without cross-run deletion or manifest
contamination, including in containers that reuse process IDs.

The examples themselves now use supported apex-centered boxes from
`default_apex_box()`, checked bundle status, explicit CSV lifecycle checks, and
moderate demonstration-scale grids.  `demo3d_2` writes auxiliary point clouds
transactionally and obtains strength, Mach number, downstream magnetic change,
and conservation residuals from `shock_state_direction_checked()`.

### Expected result

All three binaries compile without warnings and exit zero.  Stderr is empty;
stdout names the expected product or current ideal-MHD diagnostic; every
declared file exists and no undeclared or staging file remains.  All serialized
rows are finite and complete, bundle zones/dimensions match the example source,
the radial and point-cloud counts match their requests, and all CSV time axes
contain 865 samples from 0 through 259,200 seconds.  Passing runs clean their
temporary directories.

Run the build-and-execute gate with:

```sh
make demo-run
```

or, after `make` has built all targets:

```sh
./output/test_swcme --test OUT07
```

`OUT07` follows `OUT01` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST07: AMPS adapter equivalence

### What is tested

PST07 verifies that `Interface1D` and `Interface3D` return the same physical
values and status as direct SWCME queries made with the identical model,
prepared state, time, position, direction, observer, and energy grid.  It checks
prepared-state model identity and configuration digest; density, proton
pressure, velocity, magnetic components and magnitude, `div(V)`, and focusing;
shock/source position, normal, compression, `theta_Bn`, fast Mach number,
normal speed, upstream density/pressure/field, source weighting, and DSA slope;
relative and SI-normalized spectra at 1, 3, 10, 100, and 1000 MeV; and the
selected observer-to-cobpoint Parker path length.

### Why it is tested

The AMPS adapter is the operational path used by particle transport.  A wrapper
that repeats a unit conversion, default, Parker derivative, connectivity
integral, or source formula can disagree with standalone SWCME even when both
components pass isolated tests.  Such drift would make transport results depend
on the call path rather than the configured physical model.

### How it is tested

The test creates matched fast-shock/source configurations and queries multiple
upstream radii in both dimensions.  Adapter fields that should be copied are
required to equal their direct production values exactly.  Pressure, focusing,
and path length are evaluated by shared production helpers; the 3-D
connectivity class and adapter both consume that common Parker implementation.
An unsupported radius must retain the direct `OUTSIDE_MODEL_DOMAIN` status.

For source equivalence, PST07 compares the adapter record field-by-field with
the direct `ShockAccelerationState` and local 3-D shock state before any CSV
formatting.  It independently applies the documented relativistic momentum
power law at a fixed energy grid and compares both relative shape and physical
SI normalization.  Finally, it performs direct and adapter observer
connectivity calls, compares root selection and path length, and verifies local
source focusing at the selected cobpoint.

### Expected result

Every copied status and value is identical, independently ordered spectrum
calculations agree within `2e-13` relative error, model/configuration provenance
is unchanged, and invalid-domain classification is preserved.  The adapter
performs record assembly only: no independent pressure, Parker path, focusing,
shock, or source-spectrum model is allowed.

Run the gate with:

```sh
make -j
./output/test_swcme --test PST07
```

PST07 follows OUT07 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## PST05: prepared-state lifetime contract

### What is tested

PST05 tests the lifetime boundary between a prepared 1-D/3-D `StepState`, its
logical `Model` owner, and the AMPS-facing `Interface1D`/`Interface3D` adapters.
It covers model destruction, move construction, move assignment, rejected move
assignment into a prepared destination, standard-container relocation, moved-
from rejection, copied state records, and asynchronous read-only evaluation.
Both direct field APIs and adapter background APIs participate in the matrix.

### Why it is tested

A cached state that secretly refers to model-owned memory can become a dangling
record when its model leaves scope or is relocated by a return value or
container. Conversely, treating a move as an independent copy can invalidate
otherwise safe states during ordinary C++ object relocation. Either defect can
produce use-after-free, nondeterministic physics, or a plausible result from the
wrong owner. PST05 converts these implicit C++ lifetime assumptions into a
public, executable ownership rule.

### How it is tested

The test first verifies that both state types retain copy/move value semantics
and that both current parameter bundles give their models non-throwing move
construction. A prepared 1-D model is move-constructed and a prepared 3-D
model is move-assigned into an unprepared target; the destination must inherit
the original identity, frozen phase, and bitwise-identical field result. Calls
through each moved-from model must return `STATE_MODEL_MISMATCH`, report both
identities, and preserve sentinel outputs.

A second assignment attempts to overwrite an already prepared destination.
The required `std::logic_error` must occur before either identity changes, and
both pre-existing states must still evaluate successfully. Separate scoped
owners are then destroyed. Their retained states must remain safely copyable
with unchanged integrity seals, while new equal direct/adapter owners must
reject them transactionally rather than resurrecting permission from matching
parameters.

For relocation, vectors are deliberately reserved for one element and grown to
two after the first model/adapter has prepared a state. The relocated first
element must retain its identity and exact direct or AMPS output. Finally,
copies of the prepared states are captured by value in `std::async` tasks while
their relocated owners stay alive; the returned fields must exactly match the
serial pre-move records. A dedicated build compiles the full production and
validation source set under AddressSanitizer:

```sh
make -j
./output/test_swcme --test PST05
make pst05-sanitize
```

LeakSanitizer is disabled in that target by default because it cannot inspect
threads under ptrace-based CI/container supervisors. An untraced host can add
leak scanning with
`make PST05_ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 pst05-sanitize`.

### Expected result

Prepared states remain valid inert values after owner destruction, with no
invalid memory access or leak. A moved-to owner accepts all of its pre-move
states and produces identical output; the moved-from object and any newly
constructed equal owner reject them deterministically. Frozen-destination
assignment changes neither object. Container relocation and asynchronous
direct/AMPS evaluations pass with exact results, and AddressSanitizer reports no
use-after-free or invalid access. An optional LeakSanitizer run reports no leak.

PST05 follows PST07 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## PST08: state ownership performance

### What is tested

PST08 measures the runtime and memory overhead of model-identity,
configuration-digest, and prepared-state-integrity validation. It separately
times `prepare_step()`, the complete state guard, direct scalar queries, AMPS
scalar background queries, a 16,384-point 1-D fast batch, and a 256-point 3-D
fast batch. Every timed category reports a median and nearest-rank p95 after
warm-up. A C++17 allocation probe verifies the warmed validation, direct, and
adapter hot paths.

### Why it is tested

PST02, PST03, and PST06 require strong validation before physics or output
mutation. Recomputing those digests in a nested evaluator—or once for every
sample in a batch—preserves correctness but can make particle coupling
unnecessarily expensive. Before PST08, the 1-D AMPS scalar path authenticated
the same state three times, and a 3-D batch re-entered checked directional shock
and geometry methods for every point. The test prevents that topology from
returning and ensures safety checks do not introduce hidden heap traffic.

### How it is tested

The ordinary registered test uses the project's optimized build. The dedicated
`pst08-performance` target rebuilds only `swcme3d.cpp`, the benchmark, and its
minimal runner with pinned C++17, `-O3`, `-DNDEBUG`, warning, and pthread flags.
Preparation uses a fresh model in each timing iteration and is never included
in per-query throughput. Cheap scalar operations execute in inner loops so
clock resolution does not dominate; representative batches execute from
preallocated coordinate and destination arrays.

The corrected paths are compared with an emulated pre-PST08 topology in the
same process. The emulation calls the exact current state validator at each
redundant outer/nested/per-sample location removed by PST08, then calls the
corrected physics path. It therefore isolates validation topology without
keeping obsolete physics, weakening state checks, or exposing a public bypass.
All return codes and numerical outputs feed a volatile checksum to prevent
dead-code elimination.

For memory, the test executable instruments scalar, array, sized, and aligned
`operator new`. Counting is enabled only after benchmark storage and iostreams
are warm and only around public production calls. This makes any new allocation
below a digest, status, model evaluator, or AMPS record assembly visible while
leaving the production library's allocator unchanged.

Run both forms with:

```sh
make -j
./output/test_swcme --test PST08
make pst08-performance
```

### Expected result

Complete state-validation p95 is no more than 10 microseconds. Direct/AMPS
scalar p95 budgets are 20/25 microseconds in 1-D and 500/600 microseconds in
3-D. Corrected p95 batch targets are 5 milliseconds for 16,384 1-D points and
100 milliseconds for 256 3-D points; one validation contributes at most 5% of
either batch runtime. The corrected 1-D direct median must be at least 20%
faster than the emulated two-validation path, and the AMPS median at least 40%
faster than the emulated three-validation path. Warmed validation and direct/
AMPS queries allocate exactly zero bytes through C++ allocation APIs.

PST08 follows PST05 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## OUT08: strict warning writer build

### What is tested

OUT08 validates compilation of the shared status/output layer, the header-only
1-D writer, the compiled 3-D writers, and `demo1d`, `demo3d_1`, and `demo3d_2`
under the project's strict warning policy. It builds independent debug
(`-O0 -g3`) and optimized (`-O3 -DNDEBUG`) products. The policy enables extra,
pedantic, level-two format, format-security, conversion, sign-conversion, and
shadow diagnostics and promotes every selected warning to an error.

The registered runtime half also exercises the two writer entry points. A
literal containing `100%`, `%s`, and `%zu` must be copied byte-for-byte without
format interpretation; a separate typed formatted record must render the
expected integer and `std::size_t` values and complete the full checked stream
lifecycle.

### Why it is tested

Variadic output code sits at a hazardous boundary: a format/argument mismatch
is undefined behavior, runtime text used as a format creates a format-string
vulnerability, and implicit narrowing can corrupt large mesh or grid counts.
Optimized compilation can also reveal output values that are only conditionally
assigned. Ordinary `-Wall` builds did not inspect calls through the project's
custom formatter and previously diagnosed the 3-D volume/face values only in
some optimized compiler configurations. OUT08 makes these defects deterministic
CI failures.

### How it is tested

`CheckedTextFile::print()` is annotated as a GCC/Clang printf-like function, so
the compiler checks its format string and variadic arguments exactly as it
checks `printf`. `write_literal()` is a compile-time-sized array overload for
records requiring no substitution; it bypasses formatting and shares the same
raw checked-write implementation. The in-memory OUT08 probe proves both paths
produce the intended bytes without relying on production formatter constants.

The 3-D bundle and standalone-face writers initialize every scalar and call the
checked field/divergence interfaces. If evaluation fails, the writer attaches
the global row index, cancels the staging transaction, and returns the original
physics status before formatting the row. Mesh connectivity is checked while
signed, then converted once to `std::size_t`; array-axis and demonstration
index-to-time conversions are explicit. No pragma or command-line suppression
is used.

The phony target always compiles fresh isolated products, preventing an object
created with normal flags from satisfying the gate accidentally:

```sh
make -j
./output/test_swcme --test OUT08
make out08-strict CXX=g++
# In the supported Clang CI job:
make out08-strict CXX=clang++
```

The local compiler is selected through the standard overridable `CXX`
variable. CI runs the same target separately for every supported compiler; a
host does not need all compiler families installed for one invocation.

### Expected result

Both debug and optimized builds complete with zero warnings, all six
demonstration binaries/objects are produced in their isolated directories, and
the debug and optimized runtime probes each pass three assertions. Any format,
format-security, implicit conversion, sign conversion, shadowing, extra, or
pedantic diagnostic fails `make`. The literal record is preserved exactly, the
formatted record matches its expected bytes, and open/close state is complete.
No suppression is accepted unless a future exception is narrowly scoped,
documented at its call site, and covered by this test.

OUT08 follows PST08 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include its
registered runtime assertion through `@ALL`. CI and release validation must run
`make out08-strict` in addition to the selected runtime profile because a
running executable cannot verify the flags used to compile itself.

## Python campaign manager and reproducible run artifacts

`run_tests.py` is the manager for multi-test development gates and validation
campaigns.  It does not contain independent SWCME physics: each C++ validation
is still executed through the single `output/test_swcme` registry executable.
Builds also produce `output/sep_reference`, a public-interface consumer used for
reference histories and sweep probes.

Typical commands from `swcme/test` are:

```sh
python3 run_tests.py --list
python3 run_tests.py --profile SMOKE
python3 run_tests.py --profile ROUTINE
python3 run_tests.py --all
python3 run_tests.py --test SEP03
python3 run_tests.py --profile EVENT --event-config event_config.example.json
python3 python/test_run_tests.py
```

The named profiles are version-controlled text files in `profiles/`:

- `SMOKE` is the short development gate, including prepared-state safety,
  checked output failures, transactional commit, model-domain preflight,
  BoxSpec/mesh validation, independent parsing, executable demonstrations,
  state-ownership performance, strict writer record semantics, and smoothing-
  width rejection/preservation;
- `ROUTINE` is the broad deterministic gate and excludes `MSH05` and `CON05`;
- `FULL` expands to every test in the C++ registry;
- `EVENT` first runs FULL and then executes the supplied event-analysis JSON.

`--no-build` is available when the executables are already current;
`--build-only` compiles the test and reference tools without running a campaign.
`--reference-export` forces the default SEP history export and
`--no-reference-export` disables it. FULL and EVENT export it by default.

Each campaign creates an isolated output directory containing:

```text
manifest.json
summary.json
summary.csv
logs/<TEST_ID>.log
sep_reference.csv              # when reference export is enabled
sweeps/*.csv + per-case logs   # when configured
convergence/*.json             # when configured
comparisons/*.json             # when configured
plots/*                        # optional
```

The manifest records the random seed, selected tests, git metadata, compiler
path/version/flags, host and Python details, visible MPI/OpenMP environment,
source-tree SHA-256, complete resolved SWCME+SEP model configuration and hash,
and the event JSON/hash.  Use `--resolved-config FILE` when a driver has already
written the exact resolved parameter block; otherwise the manager obtains the
default block from `sep_reference --print-manifest`.

Campaign exit codes are part of the automation contract:

```text
0  all required tests/analyses passed
1  validation, reference export, or EVENT analysis failed
2  command-line/configuration error
3  build failure
```

### EVENT JSON

`event_config.example.json` is executable and demonstrates a parameter sweep.
The event object may contain `sweeps`, `convergence`, `comparisons`, and `plots`.
A sweep defines a Cartesian product of parameter arrays and a command template.
Any command token can use a sweep field plus `{seed}`, `{root}`, `{test_dir}` or
`{output_dir}`.  `metric_regex` captures a scalar from stdout; optional
`metric_min`/`metric_max` turn the captured quantity into an acceptance gate.

A convergence entry either supplies explicit positive `x` and `error` arrays or
references a previous sweep through `sweep`, `x_parameter`, and `error_metric`.
The manager fits the slope of `log(error)` versus `log(x)` and checks optional
`min_order`/`max_order` limits.

A comparison entry names model/reference CSV files, optional key/key tolerance,
and one or more columns with absolute and/or relative tolerances.  A plot entry
selects a CSV, one x column and one or more y columns.  Plotting is optional;
when Matplotlib is unavailable a non-required plot is reported as SKIP rather
than invalidating the physics campaign.

The Python manager itself is regression-tested in `python/test_run_tests.py`.
Those tests cover profile expansion against the live C++ registry, exact
second-order convergence fitting, Cartesian sweeps/metric capture, CSV
comparison and EVENT aggregation, manifest/report creation, and deterministic
CLI configuration-error handling.

## Validation classifications

- `COMMON`: contracts shared by both models, with both public paths exercised.
- `1D`: behavior specific to `swcme1d`.
- `3D`: behavior specific to `swcme3d`.
- `1D<->3D`: direct equivalence or consistency between the two implementations.

The current registry uses all four classifications.  In particular, the
`1D3D01`-`1D3D03` and `SEP03` gates use direct dimensional-equivalence fixtures
where appropriate.

## CFG01: centralized configuration rejection and physical-range validation

`CFG01` is a `COMMON` test of the production configuration boundary.  Both
model classes expose a side-effect-free `validate()` method returning
`swcme::config::ValidationResult`.  Every issue contains the public field name,
a stable error code, and the violated rule.  `prepare_step()` invokes the same
validator and throws `std::invalid_argument` before any physics if the result is
invalid, so callers cannot bypass validation accidentally.

The deterministic fixtures cover zero solar-wind speed, negative density,
negative DBM Gamma, invalid Parker normalization latitude, NaN temperature,
negative smoothing widths, malformed DATA_DRIVEN tables, zero CME direction,
zero solar-rotation axis, invalid finite-SSE half widths, non-positive
ellipsoid axis ratios, and negative solar-rotation rate.  A multiple-error
fixture confirms that one validation call reports all independent bad fields
instead of stopping after the first.  For representative invalid cases CFG01
also calls `prepare_step()` and requires rejection.

Validation is intentionally separate from unit conversion.  The common unit
layer converts `0 km/s -> 0 m/s` exactly; CFG01 separately rejects `V_sw<=0`
because the Parker/DBM baseline requires a positive ambient wind.  This is a
permanent regression guard against the former 1-D `max(1,V_sw*1000)` behavior.
The compact shared-constant baseline remains in CFG01 to guard the immutable
AU, nominal solar radius, proton mass, vacuum permeability, and Boltzmann
constant used by both interfaces.

## CFG02: centralized unit-conversion and dimensional-consistency test

`CFG02` validates `swcme_units.hpp`, the single production source for unit
conversion.  Forward and round-trip checks cover km/s, nT, cm^-3, km^-1, AU,
nominal solar radii, hours, and degrees.  Exact decimal/defined conversions are
checked exactly where possible; short floating-point chains use
`64 * std::numeric_limits<double>::epsilon()`.  No physics tolerance is used.

The test then constructs valid 1-D and 3-D models and confirms that their
prepared ambient speed, one-AU density, and one-AU magnetic-field magnitude
agree to roundoff.  An independent SI calculation of
`V_A=B/sqrt(mu0*rho)` is compared with the shared production Alfvén-speed
helper for both model paths.  Zero wind speed is deliberately *not* passed to
`prepare_step()` in CFG02: it is a valid conversion input but an invalid model
configuration and is therefore covered by CFG01.

The model wrappers now use the same conversion helpers when building their SI
state, so CFG02 no longer audits duplicated `1e3`, `1e6`, `1e-9`, or `1e-3`
conversion literals in 1-D versus 3-D.  A failure should be investigated as a
unit-contract or wrapper-integration defect; reference values or tolerances
must never be changed merely to obtain PASS.

## CFG03: smoothing width policy

### What is tested

CFG03 tests all three public region-smoothing widths in both 1-D and 3-D at
zero, representative interior values, their simultaneous exact 90-percent
limits, one ULP above each limit, and non-finite values. It also submits all
three oversized widths together, checks an oversized width while it is dormant
under `SHOCK_ONLY/SOURCE`, and audits effective widths at the 1-D shock, 3-D
apex, and an arbitrary local 3-D radius. At the exact-limit configuration it
checks that neighboring transition intervals remain separated and each
smoothstep center has a finite blend of one half.

### Why it is tested

The old region constructor accepted any non-negative width and silently used
`min(requested, 0.90*layer)` during every boundary calculation. A campaign
manifest therefore recorded one configuration while the field evaluator used
another, and the amount of hidden adjustment could depend on the local 3-D
shock radius. This is scientifically unsafe because smoothing directly controls
resolved compression and `div(V)`, hence SEP acceleration and adiabatic energy
change.

### How it is tested

The shared policy computes the maximum self-similar fractions from sheath
fraction `f_s` and ejecta fraction `f_e`:

```text
w_shock,max = 0.90 f_s
w_LE,max    = 0.90 min(f_s, f_e)
w_TE,max    = 0.90 f_e
```

Central validation compares every stored width with its limit, regardless of
whether the selected mode currently uses it. Equality must validate and prepare
successfully. `std::nextafter(limit,+infinity)` must produce a field-specific
`OUT_OF_RANGE` issue and make `prepare_step()` throw before constructing any
state. NaN and positive/negative infinity must remain `NON_FINITE`, not overlap
errors. A simultaneous invalid request must return all three issues in one
validation result.

For valid configurations, the test independently recomputes
`width_fraction*local_radius` and compares it with each prepared boundary at a
small roundoff-only tolerance. This would fail decisively if the historical
90-percent `min()` returned. Direct calls to `locate()` verify the most
restrictive accepted geometry retains finite, ordered transition regions and
canonical center weights.

Run CFG03 with:

```sh
make -j
./output/test_swcme --test CFG03
```

### Expected result

Zero, ordinary, and exact-limit inputs are accepted without modification in
both models. Each one-ULP excess is rejected before preparation with its exact
field name; non-finite classifications are preserved; all fields are reported
for a combined conflict; and dormant invalid widths are rejected. Every valid
effective width equals the recorded requested fraction times the local shock
radius to floating-point roundoff. Transition windows remain non-overlapping,
ordered, finite, and normalized.

CFG03 follows OUT08 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## CFG04: configured radius domain

### What is tested

CFG04 tests the DBM/ballistic reference radius, every data-driven PCHIP radius,
the 3-D connectivity source/search radius, and observer radii against the
inclusive `1.05 R_sun` analytical-domain boundary. Each input is exercised one
ULP below, exactly at, and one ULP above the boundary. The test also covers a
non-finite reference radius, sub-domain knots in every vector position, a
mixed-validity table, repeated times, and equal or reversed source-observer
ordering.

### Why it is tested

Previously `r0_Rs` and data-driven knots were required to be merely positive.
Such a configuration could validate and then fail only after interpolation or
shock preparation entered the unsupported Parker/Leblanc domain. Connectivity
similarly accepted a positive source radius below the model boundary and could
convert an invalid observer or empty search interval into an ordinary
`Disconnected` result. Those late or ambiguous failures obscure the actual
configuration defect.

### How it is tested

The test uses `std::nextafter` around the shared `MIN_RADIUS_RS` and
`MIN_RADIUS_M` constants, avoiding arbitrary decimal offsets. It invokes both
side-effect-free model validators and `prepare_step()` in 1-D and 3-D. Each bad
PCHIP radius must appear as `data_radius_Rs[index]`; mixed-invalid input must
report all primary knot defects in one pass while retaining independent table
ordering diagnostics. Accepted boundary configurations are prepared at the
first knot and their cached radius is compared exactly with the SI boundary.

For connectivity, the test calls the production observer solver with source
radii below, at, and above the boundary and then with equal/reversed bounds. It
also checks both observer-scope APIs and the connectivity status at the three
observer boundary points.

Run CFG04 with:

```sh
make -j
./output/test_swcme --test CFG04
```

### Expected result

Every sub-domain model radius is rejected before preparation with its public
field or knot index. Exact and above-boundary model configurations prepare
without clipping. A sub-domain connectivity source radius and reversed
source-observer ordering return `InvalidConfiguration`; equal radii form a
valid one-point closed interval. A sub-domain observer returns
`InvalidObserver`; valid boundary inputs are never mislabeled as a domain
failure. CFG04 follows CFG03 in `SMOKE` and is included in every `@ALL` profile.

## DEF01-DEF04: canonical defaults, scope, and resolved metadata

`swcme_defaults.hpp` owns every dimensionality-independent default used by the
1-D and 3-D public parameter structures.  These tests are release guards: a
change to a science default must be deliberate, documented, and changed in one
place rather than drifting independently between the two interfaces.

- `DEF01` compares all shared default ambient, kinematic, source, region,
  smoothing, sheath, and ejecta values between default-constructed 1-D and 3-D
  `Params`.  It also checks that the DBM reference radius and 1-AU density come
  from the canonical default constants.
- `DEF02` verifies the declared science convention: `SHOCK_ONLY + SOURCE` maps
  to `CONTROLLED_SEP_PRE_SHOCK`, the default 3-D geometry is finite `SSE`, the
  Parker normalization is total positive |B| at the equatorial 1-AU reference,
  radial polarity is outward, and the baseline frame name is stable.
- `DEF03` verifies observer-local scope gating.  A pre-shock observer is in the
  controlled Parker scope; after the modeled front reaches that observer the
  status is explicitly out of scope.  A late-time observer outside a finite SSE
  angular cap remains in scope because no front exists on that ray.
- `DEF04` verifies deterministic complete resolved-configuration serialization.
  Common fields, model-specific geometry, data-table counts, compatibility
  parameters, derived scope, and explicit event overrides must all appear.  The
  serialization is repeated byte-for-byte to guard deterministic campaign
  manifests/hashes.

Tests that exercise the optional `FULL_ICME/RESOLVED_COMPRESSION` path set both
options explicitly; they no longer inherit that behavior from defaults.  This is
important because a verification fixture should declare when it is testing an
optional phenomenological model rather than the controlled SEP baseline.

## DEN01: Leblanc density normalization at the reference distance

DEN01 is a `COMMON` physics/implementation test of the defining normalization
property of the production Leblanc, Dulk & Bougeret (1998) density profile. It
does not merely compare the two models with each other: an independent scalar
reference evaluates

```text
n_L(r) = A r^-2 + B r^-4 + C r^-6,
A = 3.3e5, B = 4.1e6, C = 8.0e7,
```

where `r` is heliocentric radius in nominal solar radii and the result is in
`cm^-3`. The test calculates `r_ref_Rs = AU_m / Rs_m` from the adopted constants
rather than maintaining a rounded AU-to-solar-radius constant. For each model,
the intended normalization is `n(r) = S n_L(r)`, with
`S = n_ref / n_L(1 AU)` and a fixed production reference radius of one AU.

The positive reference-density fixtures are 1, 5, 8, 20, and 100 `cm^-3`.
Each fixture constructs a fresh production model and step cache, evaluates the
public 1-D radial or 3-D Cartesian field interface at one AU, and checks

```text
abs(n_model(1 AU) / n_ref - 1) < 1e-12.
```

The acceptance criterion is strict and must not be weakened to turn a mismatch
into PASS. Production evaluators return density in `m^-3`; the test reports in
`cm^-3` for direct comparison with the configured normalization while retaining
the actual SI evaluation path.

The 3-D coverage evaluates `(+AU,0,0)`, `(0,+AU,0)`, `(0,0,+AU)`, and a
normalized non-axis-aligned `(1,2,3)` direction to verify that the ambient
profile is radial. For every reference density, DEN01 also compares 1-D against
3-D density and inferred scale factors at roundoff level. The same reference
location is supplied both as the production AU and as `(AU_m/Rs_m)*Rs_m` to
check AU/solar-radius coordinate equivalence without repeating CFG02.

Both production implementations expose their prepared SI coefficients `C2`,
`C4`, and `C6`. DEN01 independently recovers the nominal `A`, `B`, and `C` from
those caches and verifies the coefficients and inferred normalization factor at
`64 * double epsilon`. That roundoff tolerance covers only the short
floating-point operation chain; it is not a relaxed density-physics tolerance.
The model APIs do not expose a separately named scale factor, so inference from
`C2` is the closest production diagnostic.

The Leblanc coefficient literals and normalization arithmetic now live only in
production `swcme_solarwind.hpp`.  Both dimensional wrappers store the same
`StepState::common.solar_wind` cache and merely mirror `C2/C4/C6` for backward
source compatibility.  DEN01 still independently reconstructs the published
coefficients from those mirrors and executes an A=5, B=20, C=1 `cm^-3`
construction/re-evaluation sequence. Existing model A and B must retain their
original results after later instances are constructed, guarding against any
future cross-instance contamination in parallel AMPS use.

The production APIs currently fix density normalization at one AU, so tests of
a configurable 0.5- or 2-AU normalization radius are reported as `SKIP` rather
than supported by a fictitious test-only API. Values at 0.5 and 2 AU are checked
only for finite positivity as optional diagnostics. Detailed verification of
the radial `r^-2 + r^-4 + r^-6` behavior belongs to DEN02, not DEN01.

DEN01 prints the reference radius, each nominal analytical term, unscaled
density, expected and inferred scale, production density, absolute error,
normalization residual, tolerance, direction coordinates, and every individual
PASS/FAIL. Its structured text metrics include fixture/pass/fail/skip counts,
maximum normalization residual, maximum 1-D/3-D difference, coefficient
mutation status, and cross-instance contamination status.

A DEN01 failure may indicate incorrect radial units or coefficient powers, an
incorrect normalization factor, a `cm^-3`/`m^-3` conversion error, inconsistent
1-D and 3-D implementations, mutable global/static normalization state, or
accidentally modified nominal coefficients. The independent analytical helper
must remain separate from the production routines under test; references and
tolerances must never be changed solely to make a production failure pass.

## CFG01: centralized configuration rejection and physical-range validation

`CFG01` now exercises the production validation API rather than reporting the
configuration layer as a skip.  Both model classes expose a side-effect-free
`validate()` method returning `swcme::config::ValidationResult`.  Every issue
contains the public field name, an error code, and the violated rule.
`prepare_step()` calls this same validator and throws `std::invalid_argument`
before any physics if the result is invalid.

The deterministic fixtures cover zero solar-wind speed, negative density,
negative DBM Gamma, invalid Parker normalization latitude, NaN values, negative
smoothing widths, malformed DATA_DRIVEN tables, zero CME direction, zero solar
rotation axis, invalid finite-SSE half widths, non-positive ellipsoid axis
ratios, and negative solar rotation rate.  A multiple-error fixture confirms
that one validation call reports all independent bad fields rather than only
the first.  For representative invalid cases CFG01 also calls `prepare_step()`
and requires rejection, proving that no caller can bypass the validation layer
by ignoring `validate()`.

Validation is intentionally distinct from conversion.  For example, the unit
layer maps `0 km/s -> 0 m/s` exactly; CFG01 separately rejects `V_sw<=0` because
the baseline Parker/DBM model requires positive ambient wind speed.  This is a
permanent regression guard against the former 1-D `max(1,V_sw*1000)` behavior.

## CFG02: centralized unit-conversion and dimensional-consistency test

`CFG02` validates the production helpers in `swcme_units.hpp` for forward and
round-trip conversions of velocity, magnetic field, number density, inverse
length, AU, solar radii, hours, and angles.  Exact decimal/defined conversions
are checked at exact or roundoff precision; no physics tolerance is used.

The test then constructs valid 1-D and 3-D models and confirms that their
prepared solar-wind speed, 1-AU density, and 1-AU magnetic-field magnitude agree
to roundoff.  Finally, an independent SI calculation of
`V_A=B/sqrt(mu0*rho)` is compared with the shared production Alfven-speed helper
for both model paths.  The previous CFG02 failure at zero wind speed is removed
by design: zero is a valid unit-conversion input and an invalid model
configuration, so it is tested in the correct layer rather than by forcing a
model setup with inadmissible parameters.

## Current SWCME unit contract

External `Params` values use:

- velocity in km/s;
- magnetic-field magnitude in nT;
- number density in cm^-3;
- launch radius in nominal solar radii;
- drag coefficient in km^-1;
- sheath/ejecta dimensions in AU normalized at 1 AU;
- cone angle in radians, despite the 3-D default being initialized from 40
  degrees; and
- temperature in kelvin.

Model evaluators and prepared states use SI meters, seconds, m/s, tesla, m^-3,
and radians. SWCME adopts `1 AU = 149597870700 m`, nominal solar radius
`6.957e8 m`, CODATA 2022 proton mass `1.67262192595e-27 kg`, vacuum
permeability `1.25663706127e-6 N/A^2`, exact Boltzmann constant
`1.380649e-23 J/K`, and the model solar-rotation convention `2.86533e-6 rad/s`.

Unit conversion is now centralized in production `swcme_units.hpp`; the 1-D
and 3-D `prepare_step()` paths no longer maintain independent km/s, cm^-3, nT,
or km^-1 conversion factors.  Conversion and physical validation are separate:
zero values remain zero under conversion, while `swcme_config.hpp` decides
whether a model parameter is admissible.  CFG02 verifies the common conversion
helpers directly and then checks that both dimensional wrappers produce the
same prepared SI state for valid physical inputs.

## PAR01-PAR03: corrected 3-D Parker magnetic field

The production 3-D Parker field now uses an explicit solar-rotation axis and a
local spherical basis at every evaluation point.  The azimuthal direction is

```text
e_phi = (Omega_hat x e_r) / |Omega_hat x e_r|
```

and the local pitch is proportional to
`|Omega_hat x e_r| = sin(theta_local)`.  This corrects the previous
implementation, which used `(Omega_hat x e_r) x e_r` (a meridional direction)
and one global `sin_theta` for the entire 3-D domain.

`Params::solar_rotation_axis` specifies the global solar-rotation axis and is
normalized when a `StepState` is prepared.  A zero or non-finite axis is
rejected because it cannot define a Parker azimuthal direction.  The existing
`Params::sin_theta` field is retained only for backward-compatible
normalization of `B1AU_nT`: it specifies the reference `sin(colatitude)` at
which `B1AU_nT` is interpreted as total field magnitude at 1 AU.  It no longer
sets the local 3-D Parker winding.  The step cache stores `k_AU = Omega*AU/Vsw`;
each point multiplies this by its own `r_AU*sin(theta_local)`.

The pole is handled analytically.  When the radial direction is parallel to
the rotation axis, `sin(theta_local)=0`, so `B_phi=0` and the field is purely
radial.  The implementation does not manufacture an arbitrary azimuthal unit
vector at the coordinate singularity.

The following deterministic 3-D tests were added:

- `PAR01` — equatorial Parker-vector orientation.  With the rotation axis along
  +Z and the point at +X, it requires the outward radial component to lie along
  +X, the Parker azimuthal component to lie along -Y for the current polarity
  convention, and the meridional Z component to vanish.  It also verifies the
  requested total-field normalization at the reference latitude.
- `PAR02` — arbitrary latitude and arbitrary rotation axis.  It compares the
  production vector with an independent analytical Parker reference away from
  the equator and repeats the check after rotating the solar axis, preventing a
  hard-coded +Z implementation from passing accidentally.
- `PAR03` — polar-limit regularity.  It verifies the exact north and south
  rotation poles and a near-pole point, requiring the transverse field to tend
  continuously to zero without NaN, Inf, or an arbitrary transverse direction.

The analytical reference in `3d/test_parker.cpp` is intentionally independent
of the production Parker helper.  Component comparisons use roundoff-level
absolute tolerances scaled to the expected Tesla magnitude.  These tests do not
cover numerical solenoidality or field-line tangent/path-length validation;
those remain separate planned tests (`PAR04` and later) so the individual
validation requirements remain diagnostically focused.

## GEO01-GEO08: corrected finite shock geometry, normals, and flank speed

The 3-D shock geometry has been reworked so that the production model no longer
uses the former `ConeSSE` approximation

```text
R(theta) = R_apex cos(theta)^m
```

with a radial normal and an artificial clamped radius beyond the configured
half width.  `ShockShape::SSE` now denotes a true finite self-similar-expansion
spherical cap.  `ShockShape::ConeSSE` is retained only as a source-compatible
enum alias and has the same corrected SSE semantics; `flank_slowdown_m` remains
in `Params` only for source/input compatibility and is ignored by the SSE
geometry.

For an apex distance `R_apex` and angular half width `lambda`, the generating
sphere has center distance and radius

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda).
```

A heliocentric ray separated from the CME axis by `alpha` intersects the
outward front at

```text
R(alpha) = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2),
```

provided `alpha <= lambda`.  At `alpha=lambda` the ray is tangent to the
generating sphere.  Directions beyond the half width return `exists=false`;
the public geometry API also returns zero radius/normal sentinels so a caller
cannot mistake an absent surface for a physical flank.

The exact outward normal is the normalized level-set gradient

```text
n_hat = (R e_r - c e_CME) / a.
```

All supported shapes are treated as self-similar.  For a fixed surface
direction,

```text
dR/dt = V_apex R/R_apex,
V_sh,n = V_apex (R/R_apex) (e_r dot n_hat).
```

This replaces the former independent `cos(theta)^m` flank-speed factor and
also corrects the ellipsoid, whose flanks previously inherited the full apex
speed.  For a sphere the formula reduces exactly to `V_sh,n=V_apex`; at the
mathematical SSE tangent boundary the normal projection tends to zero.

`Model::shape_radius_normal()` and `Model::diagnose_direction()` now return a
boolean surface-existence flag.  Existing callers that ignore the return value
remain source-compatible, but finite-width-aware code should always test it.
The Cartesian plasma/field evaluators do so internally: outside the SSE angular
support they return the undisturbed ambient solar wind/Parker field and do not
construct a sheath, ejecta, shock normal, or magnetic amplification from a
fabricated flank.

The following deterministic tests validate the corrected geometry:

- `GEO01` — Sun-centered spherical reference: radius is direction independent
  and the outward normal equals the radial unit vector.
- `GEO02` — true SSE apex: the directional radius equals the configured apex
  distance and the normal equals the CME propagation direction.
- `GEO03` — SSE tangent flank: the boundary point at `alpha=lambda` remains a
  finite valid surface point; the ray is tangent and the normal radial
  projection tends to zero.
- `GEO04` — finite-width enforcement: any direction beyond the half width
  returns `exists=false`, zero geometry sentinels, `rc=1`, and `V_sh,n=0` from
  the diagnostic API.
- `GEO05` — SSE level-set residual: returned surface points satisfy
  `|x-c e_CME|=a` to near machine precision over multiple angles/azimuths.
- `GEO06` — analytical normal validation: SSE normals are compared with an
  independent finite-difference gradient of the dimensionless spherical level
  set; the ellipsoid normal is independently checked against its level-set
  gradient.
- `GEO07` — normal-speed validation: the reported `V_sh,n` for sphere, SSE, and
  ellipsoid is compared with the normal projection of centered finite-
  difference surface motion at neighboring times.
- `GEO08` — rotational covariance: rotating the CME axis, solar axis, and query
  direction together must rotate the normal while leaving radius, normal
  speed, and scalar physical compression unchanged.

The analytical/reference calculations in `3d/test_geometry.cpp` intentionally
do not call the production geometry helper.  They are independent checks of
surface radius, level-set membership, normal orientation, and motion.  Test
tolerances are roundoff- or finite-difference-level tolerances appropriate to
each calculation and must not be loosened simply to make a production result
pass.

The geometry tests remain diagnostically focused on shape, normals, finite
angular support, and normal speed.  Shock existence/compression/downstream
physics is covered independently by `SHK01`-`SHK12` below, while the now-
corrected triangulation/topology is qualified independently by `MSH01`-`MSH05`.

## MSH01-MSH05: shock-surface mesh topology, quality, and area sampling

The production mesh is no longer a rectangular theta-phi array.  It stores one
apex, unique periodic rings, and either one finite-SSE boundary ring or one rear
pole for a closed Sphere/Ellipsoid.  The seam is represented by wrapped
connectivity, so there is no second copy of the `phi=0` vertex at `phi=2*pi`.
Triangle metrics reject repeated indices, scale-aware degenerate area, and
inward winding as `INVALID_MESH`; validation must never make a bad mesh pass by
filtering cells after construction.

- `MSH01` — **nondegeneracy**.  Builds minimum, typical, and fine Sphere/SSE
  meshes, including narrow and broad SSE caps.  Every triangle must have finite
  positive area and `compute_triangle_metrics()` must accept the complete mesh.
  The test prints minimum/median/maximum area for each fixture.
- `MSH02` — **outward cell orientation**.  For Sphere, rotated Ellipsoid, and
  rotated SSE meshes, each triangle cross-product normal is compared with the
  analytical surface normal evaluated in the centroid direction.  Every dot
  product must be positive; the worst alignment is reported.
- `MSH03` — **surface-area convergence**.  Refines `nTheta` and `nPhi` by factors
  of two and sums physical triangle area.  The Sphere reference is `4*pi*R^2`.
  For true SSE with apex radius `R` and half width `lambda`, the independent
  translated-sphere reference is `2*pi*R^2*sin(lambda)^2/(1+sin(lambda))`.
  Error must decrease monotonically, the final observed order must exceed 1.5,
  and the fine-grid error must be below 0.5%.
- `MSH04` — **unique apex / periodic seam**.  Inspects connectivity rather than
  only coordinates.  A finite SSE cap must have `1+nTheta*nPhi` vertices,
  `nPhi*(2*nTheta-1)` triangles, no duplicate coordinate pairs, exactly
  `nPhi` boundary edges, manifold edge incidence, apex valence `nPhi`, disk
  Euler characteristic one, and explicit last-to-first ring adjacency.
- `MSH05` — **area-weighted source-patch sampling**.  Verifies the production
  cumulative-area table is strictly increasing and ends exactly at one.  Three
  fixed RNG seeds each draw 200,000 cells through `sample_triangle_by_area()`.
  Aggregated counts are compared with exact physical area probabilities using a
  chi-square statistic and independent incomplete-gamma p-value calculation;
  every deterministic fixture requires `p > 1e-3`.

Run this mesh gate directly with:

```sh
./output/test_swcme --test MSH01
./output/test_swcme --test MSH02
./output/test_swcme --test MSH03
./output/test_swcme --test MSH04
./output/test_swcme --test MSH05
```

The random-number generator intentionally remains outside the production mesh
class.  Production code provides only the area CDF and deterministic mapping
from a caller-supplied uniform variate to a triangle, so AMPS controls seeds and
parallel RNG policy without being able to accidentally revert to uniform-by-cell
sampling.


## DIV01-DIV03: velocity-divergence treatment

These tests qualify the production divergence operators used by SEP adiabatic
energy change.  They enforce the rule that a radial identity may only be used
for a radial velocity field and that the general 3-D operator must demonstrate
its numerical order on an independent manufactured solution.

- `DIV01` — **analytical constant radial wind**.  Configures `SHOCK_ONLY` in
  both 1-D and 3-D, samples several radii and several unrelated Cartesian
  directions, and requires the production result to equal `2*Vsw/r` to
  roundoff.  Directional dependence is a failure.  This test also guarantees
  that the baseline calculation does not reintroduce finite-difference noise.
- `DIV02` — **manufactured nonconstant radial flow**.  Uses
  `Vr=V0(1+a x+b x^2)`, `x=r/r0`, for which the symbolic spherical divergence is
  known exactly.  `swcme::divergence::radial_terms()` must agree at better than
  `1e-10` relative error.  The same test then differentiates representative
  production shock/sheath/LE/TE velocity samples independently and confirms the
  analytical `RadialVelocityState::d_velocity_dr_s_inv` follows the actual
  transport profile.
- `DIV03` — **general 3-D Cartesian divergence**.  A cubic manufactured vector
  field is evaluated at four successively halved steps.  Because centered
  differences are not exact for the cubic terms, the observed error must show
  second-order convergence.  The test additionally forces the production 3-D
  Cartesian operator on an off-axis constant radial wind and verifies
  convergence toward the exact `2*Vsw/r` result.

Run the divergence gates directly with:

```sh
./output/test_swcme --test DIV01
./output/test_swcme --test DIV02
./output/test_swcme --test DIV03
```

A passing implementation must not obtain `FULL_ICME` divergence by projecting
`V` onto `e_r` and differentiating only along the ray.  That formula drops
non-radial Rankine-Hugoniot velocity and angular gradients of a finite shock
surface.  Conversely, `SHOCK_ONLY` should not be made noisier by forcing the
Cartesian finite-difference operator when the exact analytical result is known.



## KIN01-KIN09: shared CME/shock-apex kinematics

`swcme_kinematics.hpp` is the production source of apex radius and speed for
both `swcme1d` and `swcme3d`.  The kinematics tests are classified `COMMON` and
exercise both the common SI solver and the two public dimensional wrappers.
The independent DBM reference in `core/test_kinematics.cpp` evaluates the
closed-form equations directly and does not call the production DBM helper.

The production DBM solves

```text
DeltaV0 = V0 - Vsw
DeltaV(t) = DeltaV0 / (1 + Gamma |DeltaV0| t)
R(t) = R0 + Vsw t
       + sign(DeltaV0) log(1 + Gamma |DeltaV0| t) / Gamma.
```

This sign-aware form is required for both fast and slow CMEs.  `Gamma=0` uses
the exact ballistic solution.  For very small `x=Gamma*|DeltaV0|*t`, the
production code evaluates the logarithmic distance through a short
`log1p(x)/x` series; `KIN04` explicitly straddles that numerical branch
boundary to guard against a discontinuity.

The data-driven mode uses monotone PCHIP radius interpolation.  The input time
table must be strictly increasing and radius must be nondecreasing.  The
interpolant passes through each knot exactly and its derivative is returned as
the apex speed.  Out-of-range queries return `OUTSIDE_TIME` by default;
ballistic endpoint continuation is available only when requested explicitly.

The tests are:

- `KIN01` — fast-CME DBM versus the independent closed form at multiple times,
  including exact 1-D/3-D wrapper equivalence;
- `KIN02` — slow-CME sign-aware branch; verifies acceleration toward `Vsw`
  without clipping or overshoot and checks both wrappers;
- `KIN03` — exact `Gamma=0` ballistic radius/speed in the common, 1-D, and 3-D
  paths;
- `KIN04` — small-`Gamma` analytical accuracy and continuity across the
  series/direct-log numerical branch;
- `KIN05` — long-time stability, monotonic radius, finite state, and monotonic
  approach of speed toward the ambient wind for fast and slow CMEs;
- `KIN06` — data-driven PCHIP exactness at every height-time knot, derivative
  consistency, and 1-D/3-D use of the same trajectory;
- `KIN07` — dense-grid monotonicity, nonnegative propagation speed, and absence
  of cubic overshoot between data knots; and
- `KIN08` — explicit out-of-time status, optional ballistic continuation, and
  rejection of duplicate-time or decreasing-radius tables; and
- `KIN09` — forward/backward extrapolation domain, independent formulas,
  radial crossings, a slow-DBM turning point and denominator pole, overflow
  extremes, stationary continuation, and 1-D/3-D status propagation.

The DBM analytical comparisons use `1e-12` relative accuracy where specified
by the validation plan.  PCHIP knot radii are required to be exact to
`1e-13` relative.  Tolerances must not be relaxed merely to obtain PASS.

The 1-D and 3-D public `Params` structures expose the same kinematic mode plus
`data_time_s`, `data_radius_Rs`, and `data_extrapolation`.  The 1-D convenience
method `SetDataDrivenKinematics()` sets these fields and switches the mode to
`DataDriven`.  The wrappers convert radii from solar radii to SI only when
constructing the common configuration.

A data-driven query outside the measurement interval causes `prepare_step()`
to throw a clear runtime error unless a continuation policy was selected.  This
is intentional: silent cubic extrapolation is a modeling assumption and is
not allowed to masquerade as measured/constrained kinematics.

KIN09 adds a separate result contract for an explicitly selected continuation.
Malformed configuration or a non-finite time is `INVALID_INPUT`; refusal to
leave a measured PCHIP interval is `OUTSIDE_TIME`; and a valid formula request
that crosses the radial boundary, enters a negative-speed branch, reaches the
DBM denominator pole, or overflows is `OUTSIDE_DOMAIN`. Zero speed is accepted
for the existing stationary-front and flat-knot use cases.

### KIN09 detailed validation procedure

**What is tested.** Ballistic, fast- and slow-CME DBM, and explicit PCHIP
ballistic continuation are evaluated before and after their reference or data
epochs. Fixtures include logarithmically spaced offsets, radial crossings, a
slow-DBM speed turning point, the DBM backward-time pole, maximum finite
floating-point times, and zero-slope continuation.

**Why it is tested.** An extrapolator that reports `OK` with a non-finite or
sub-domain radius can contaminate region geometry, shock timing, and every
downstream field. Treating the same event as generic invalid input loses the
important distinction between malformed configuration and a valid model
extended beyond its supported domain.

**How it is tested.** The test recomputes ballistic and endpoint-continuation
radii in long double and evaluates the sign-aware DBM formula independently
with long-double `log1p`. Valid results are compared at roundoff-level
tolerances. Points immediately around the slow-CME turning time and at the DBM
pole verify branch classification. Extreme finite time values force
intermediate or result overflow. Finally, valid public data-driven
configurations are extrapolated beyond the radial domain through both wrappers.

**Expected result.** Every representable state with `r >= 1.05 R_sun` and
nonnegative finite speed matches the independent formula. Unsupported radial,
negative-speed, pole, and overflow states return `OUTSIDE_DOMAIN` with NaN
payload fields; they are never `OK`, clipped, or confused with `OUTSIDE_TIME`.
Both wrappers propagate the status before returning a prepared state.

Run the focused test with:

```sh
make -j
./output/test_swcme --test KIN09
```

## SHK01-SHK16: fast-shock existence, Rankine-Hugoniot validation, and shock-surface ownership

The `SHK` group validates the shared production shock solver in
`swcme_shock.hpp` and its integration into both the 1-D and 3-D models.  SHK01-
SHK12 and SHK15-SHK16 exercise the shared jump physics and are classified
`COMMON`; SHK13-SHK14 are 3-D integration guards for the rule that local shock
physics belongs to the shock surface rather than to an arbitrary query point.

The solver first evaluates the upstream fast-mode speed along the shock normal.
A geometric front is not automatically a shock.  The physical condition is

```text
U1n = Vsh,n - V1.n > c_fast,
M_fast = U1n/c_fast > 1.
```

If this condition is not met, the result is `has_shock=false`, compression is
exactly one, and downstream equals upstream.  No sheath compression floor is
allowed to override this decision.

For a physical fast shock, the solver works in a frame moving with the normal
shock speed.  For a trial compression `r=rho2/rho1`, mass conservation fixes
`u2n=u1n/r`.  Tangential electric-field and tangential momentum continuity form
a 2x2 linear system for `B2t` and `u2t`; normal momentum determines `p2`.  A
bracketed solve of total-energy-flux conservation selects the compressive
fast-shock branch.  The trivial `r=1` solution is divided out of the scalar
residual so it cannot be mistaken for the physical shock root.

The individual tests are:

- `SHK01` — fast-shock existence and no-shock threshold.  Includes permanent
  1-D and 3-D regression checks proving that `sheath_comp_floor` cannot create
  a shock when `Vsh=Vsw`.
- `SHK02` — acute `theta_Bn` calculation and invariance under `B -> -B`.
- `SHK03` — strictly parallel limit against the independent gas-dynamic normal
  shock solution.
- `SHK04` — strictly perpendicular ideal-MHD benchmark against an independent
  test-side scalar reduction.
- `SHK05` — independent 80-digit benchmark for the complete eight-equation
  oblique ideal-MHD jump system across Mach number, `theta_Bn`, plasma beta,
  gamma, density/temperature, field strength and polarity, and tangential flow.
- `SHK06` — mass-flux conservation.  It also samples the production 3-D field
  immediately downstream and verifies that density, velocity, and magnetic
  field approach the exact RH downstream state instead of the ambient wind.
- `SHK07` — continuity of the normal magnetic-field component.
- `SHK08` — tangential ideal-MHD electric-field continuity.
- `SHK09` — vector momentum-flux conservation.
- `SHK10` — total ideal-MHD energy-flux conservation.
- `SHK11` — physical admissibility: positive downstream pressure/density,
  compressive branch, strong-shock bound for gamma=5/3, and entropy increase.
- `SHK12` — two-sided near-Mach-one sweep, independent full-state continuation
  reference, continuous compression limit, root diagnostics, and explicit
  classification of sub-resolution or near-singular weak branches.
- `SHK13` — shock-state independence from arbitrary query radius.  The test
  calls the legacy scalar wrapper with query radii ranging from well inside to
  far outside the front and requires identical compression, normal speed, and
  `theta_Bn`.  It also proves the fixture is sensitive by confirming that the
  ambient density at those radii differs materially from the density at the
  physical shock surface.
- `SHK14` — canonical-state consistency across `diagnose_direction()` and shock
  mesh nodes.  Every reported radius, normal, compression, and normal speed is
  compared with a fresh `shock_state_direction()` query.  This prevents future
  diagnostic/mesh code from reintroducing a separate shock-strength calculation.
- `SHK15` — fixed-seed 100,000-case stratified shock stress campaign with
  finite-output, status-classification, and evolutionary-branch requirements.
- `SHK16` — determinant-conditioned tangential-system validation on both sides
  of the pole, including one-ULP perturbations and bracket segmentation.

The independent parallel/perpendicular references and the general-oblique
`SHK05` reference do not call the production nonlinear shock solver.
Conservation tests recompute the conserved fluxes from the returned primitive
states.  This prevents a test from passing merely by reusing the same internal
algebra that produced the result.

### SHK05 independent oblique-shock benchmark

**What is tested.** `SHK05` exercises twelve well-conditioned evolutionary
fast shocks.  The matrix spans fast Mach number 1.5--6, `theta_Bn` 15--75
degrees, plasma beta 0.1--5, gamma 1.4--5/3, upstream number density
1.5--25 cm^-3, magnetic magnitude 2--15 nT, nonzero y/z tangential flow,
non-coplanar magnetic azimuths, and both field polarities.  Pressure and density
vary independently, so the upstream temperature range is varied as well.  The
test compares fast Mach number, compression, downstream density and pressure,
all three downstream velocity components, all three downstream magnetic-field
components, entropy ratio, downstream fast Mach number, branch label, and all
five normalized production conservation residuals.

**Why it is tested.** Earlier oblique coverage derived its expectations from
the same reduced relations used by production.  Such a test can confirm
internal consistency while missing a shared algebra error or selection of a
non-evolutionary nonlinear root.  `SHK05` supplies an external numerical oracle
for the general oblique case between the independent parallel (`SHK03`) and
perpendicular (`SHK04`) limits.

**How it is tested.** The reviewed v1 values live in
`reference/shk05_oblique_v1.hpp`.  The audit-only generator
`reference/generate_shk05_oblique_v1.py` uses standard-library `Decimal`
arithmetic at 80-digit precision and solves density, pressure, three velocity
components, and three magnetic components simultaneously from the complete
Rankine-Hugoniot system.  It uses a numerically formed full 8x8 Jacobian and
independent Gaussian elimination, rather than the production compression
polynomial/reconstruction.  Six starting compressions expose competing roots.
An eight-step, ten-percent Mach-continuation round trip must return to the same
root within `1e-45`; compression, entropy, and upstream/downstream fast and
Alfven characteristics then independently identify the unique evolutionary
fast branch.  The header records convergence counts, root count, maximum
high-precision residual, and minimum Newton pivot for diagnostic retention.
Normal builds depend on but never execute the generator, so production changes
cannot silently rewrite the oracle.

**What is expected.** All twelve production solves return `Solved`; the stored
reference residual is below `1e-50`, at least two distinct starting seeds reach
the physical solution, exactly one admissible root is found, and the recorded
minimum pivot exceeds `1e-10`.  Production and reference primitive quantities
must agree within `2e-9` relative error (tighter than the required `1e-7`), the
independent and production branch classifications must both be evolutionary
fast, and every normalized conservation residual must meet its production
acceptance threshold.  A failure retains the fixture case name plus branch,
conditioning, root-count, and residual diagnostics.  `SHK05` follows `CON10`
in the priority-ordered `SMOKE` profile and is also included in all `@ALL`
profiles.

### SHK12 near-Mach-one shock limit

**What is tested.** `SHK12` covers the two-sided limit about `M_fast=1` for
plasma beta 0.01, 0.1, 1, and 10; `theta_Bn` values from 1 to 89 degrees; gamma
1.4 and 5/3; nonzero tangential flow; and non-coplanar magnetic fields.  The
frozen portion contains 48 well-conditioned states at six logarithmically
spaced Mach excesses from 0.5 through `1e-5`.  For each state the test checks
Mach excess, compression, density, pressure, all velocity components, all
magnetic components, iteration count, and final compression bracket.  A
separate 32-family runtime matrix approaches the threshold from both sides down
to `|M_fast-1|=1e-12` and includes nearly parallel, low-beta configurations
whose tangential system is close to singular.

**Why it is tested.** The former three-point SHK12 exercised only one beta,
angle, and gamma, and only above the threshold.  It could not detect a physical
weak shock being reported as `NO_SHOCK`, an intermittent bracket failure, or
selection of a finite-compression outer root.  The original outermost-bracket
policy did exhibit that last defect, returning compression approximately 3--6
for some low-beta states as `M_fast` approached one.  Weak shocks require an
explicit separation between physical classification and numerical resolution
because compression alone cannot distinguish those failures.

**How it is tested.** The reviewed fixture is
`reference/shk12_near_mach_v1.hpp`; its audit-only generator is
`reference/generate_shk12_near_mach_v1.py`.  The generator imports only the
independent Decimal conservation equations and linear algebra shared with the
SHK05 reference tooling, then solves density, pressure, three velocities, and
three magnetic components together at 80-digit precision.  Multiple initial
compressions establish the evolutionary-fast root at `M_fast-1=0.5`, and
fifteen continuation points follow that same branch to `1e-5`; six points per
family are frozen.  The production solver retains continuous nontrivial
compression brackets in ascending order, solves and independently classifies
them, publishes the accepted bracket endpoints/width and iteration count, and
remains inside a deliberately broad weak-branch continuity envelope.
Normal builds never regenerate the fixture.

**What is expected.** Subcritical and exactly critical states return
`NO_SHOCK`, `has_shock=false`, compression one, and the unchanged upstream
state.  Positive Mach excesses through the published binary64 limit of `1e-6`
return `NUMERICALLY_UNRESOLVED_WEAK_SHOCK`, `has_shock=true`, and
`solver_converged=false`; they never become `NO_SHOCK` and never contain NaN or
infinite payload values.  A near-singular scan that can reach only a
discontinuous outer root receives the same explicit unresolved status.
Well-conditioned points from `1e-5` to order unity return `SOLVED`, converge in
at most 120 iterations, contract the bracket to the declared tolerance,
approach compression one monotonically, and match every frozen downstream
primitive component within `2e-7` relative error.  Reference residuals remain
below `1e-50`.  `SHK12` follows `SHK05` in the priority-ordered `SMOKE` profile
and is included in every `@ALL` profile.

### SHK16 near-singular tangential system

**What is tested.** `SHK16` conditions the oblique ideal-MHD tangential 2x2
system to a known determinant zero inside the supported compression interval.
It samples logarithmic offsets on both sides of that pole, including offsets
inside the documented singularity threshold, and repeats the complete shock
solve at the nominal shock speed and its two adjacent binary64 values.  The
test checks the signed determinant, singular classification, public
conditioning diagnostics, final bracket endpoints, and evolutionary-fast
branch flag.

**Why it is tested.** The energy residual is rational at the tangential-system
pole.  Skipping an invalid candidate without ending the current scan segment
can connect residuals from opposite sides and manufacture a false root.  An
arbitrary bisection contraction after a singular midpoint can likewise lose
the sign-bearing endpoint and converge to an unverified branch.

**How it is tested.** The fixture independently sets the singular compression
from `rho*u1n^2/r=Bn^2/mu0` and calculates every expected relative determinant
directly from that equation.  The chosen compression is exactly a nominal
scan abscissa.  Offsets of `1e-3`, `1e-6`, `1e-9`, and `1e-13` are evaluated on
both sides; one-ULP shock-speed perturbations then exercise the public solver.
Production records the closest signed and absolute conditioning values, breaks
the bracket history at every invalid point, and returns a distinct
`NUMERICALLY_SINGULAR` status when no verified continuous bracket remains.

**What is expected.** Determinant signs agree with the independent formula,
the `1e-13` candidates are singular while wider offsets remain nonsingular,
and all diagnostic values are finite.  One-ULP perturbations retain the same
classification.  Every outcome is either a verified `SOLVED` evolutionary
fast branch whose final bracket lies wholly on one side of the pole, or an
explicit `NUMERICALLY_SINGULAR` rejection; no unclassified bracket or wrong
branch is accepted.  `SHK16` follows `SHK12` in priority and in `SMOKE`.

### SHK15 high-count random shock stress

**What is tested.** `SHK15` executes exactly 100,000 deterministic physically
valid inputs in four explicit strata: 25,000 sub-fast states, 25,000
super-fast states within the published weak-shock resolution, 49,900 resolved
random shocks, and 100 determinant-conditioned shocks.  The matrix spans
log-uniform density from 0.01 to 100 cm^-3, magnetic strength from 0.1 to
100 nT, beta from `1e-3` to 100, the resulting temperature range, arbitrary
three-dimensional normals and tangential flows, both magnetic polarities,
gamma from 1.2 to 5/3, and broad fast-Mach excesses.

**Why it is tested.** A small benchmark matrix cannot expose rare interactions
among obliquity, plasma scales, weak roots, and determinant conditioning.
High-count randomized coverage detects NaN/Inf propagation, intermittent
bracket loss, unclassified reconstruction or conservation failures, and
selection of a non-evolutionary branch while remaining exactly reproducible.

**How it is tested.** A local SplitMix64 sequence with seed
`0x53484b31355f7631` and an explicit 53-bit floating conversion avoids
implementation-dependent standard-library distributions.  All geometric
frames are constructed by stable cross products.  Sub-fast and weak-limit
strata have exact expected statuses; resolved cases may return only `SOLVED`
or a documented numerical-limit rejection.  Conditioned cases deliberately
place a determinant zero on a production scan node.  Any unexpected outcome
prints its full index, stratum, normal, primitives, gamma, requested Mach,
shock speed, and status as a reusable minimized-fixture starting point.

**What is expected.** Every scalar and primitive component is finite.  All
100,000 cases return `SOLVED`, `NO_SHOCK`,
`NUMERICALLY_UNRESOLVED_WEAK_SHOCK`, or `NUMERICALLY_SINGULAR` consistently
with their stratum.  Solved cases are evolutionary-fast.  At least 90
conditioned scans observe a singular trial, and there are zero generic
no-bracket, invalid-reconstruction, conservation, or wrong-branch outcomes.
The test follows `SHK16` in priority and in `SMOKE`.

### SHK06 independent mass-flux conservation

**What is tested.** `SHK06` independently evaluates normal mass flux for every
`SOLVED` member of the 100,000-case SHK15 campaign and all 48 frozen SHK12
weak-shock states.  It explicitly requires substantial weak-shock coverage
below compression 1.01 and high-compression coverage above four.  The existing
3-D resolved-compression integration check remains and verifies that fields
immediately downstream approach the complete Rankine-Hugoniot state.

**Why it is tested.** A stored solver residual can agree with the algebra that
constructed a downstream state while both are wrong.  Recomputing the flux
from only the serialized primitive records catches a damaged density or normal
velocity component, frame-sign mistakes, and scale-dependent cancellation.

**How it is tested.** The test reconstructs each shock-frame velocity as
`u=V-Vsh*n`, takes its normal projection in long-double arithmetic, and forms
`rho*u_n` independently on both sides.  Their absolute difference is divided
by the larger physical flux magnitude, with only the smallest representable
long-double value as a zero guard.  It does not read `JumpResult.mass_residual`.
Failures print the fixed-seed stress record or frozen fixture identifier.

**What is expected.** At least 49,000 solved random/frozen states are checked,
including at least eight weak and eight high-compression examples.  Every
independent normalized residual is finite and no greater than the production
acceptance threshold `1e-9`.  The 3-D inner shock-layer sample continues to
match the exact downstream density, velocity, and magnetic field within its
existing `2e-7` relative tolerance.  `SHK06` follows `SHK15` in priority and
in `SMOKE`.

### SHK07 independent normal magnetic-field continuity

**What is tested.** `SHK07` independently projects the upstream and downstream
magnetic fields onto every arbitrary three-dimensional normal in the solved
SHK15 population.  It additionally creates 256 proper-rotation pairs and 256
magnetic-polarity-reversal pairs, requiring unchanged solved compression and
normal-field continuity for each transformed case.

**Why it is tested.** The ideal-MHD divergence constraint requires
`B1.n=B2.n`.  A test using only an axis-aligned normal or the production
residual could miss a component-order bug, a non-unit projection error, an
orientation-dependent reconstruction, or incorrect handling of negative
magnetic polarity.

**How it is tested.** Both projections are accumulated directly from the
serialized Cartesian components in long-double arithmetic.  Their difference
is normalized by the larger normal-field magnitude, with a `1e-12` fraction of
the total field as the perpendicular-limit guard.  A cyclic Cartesian
permutation supplies a determinant-+1 rotation; all vector inputs are rotated
together.  A separate pair negates the complete upstream field.  No production
projection or stored `normal_B_residual` is used as the oracle.

**What is expected.** At least 49,000 solved cases pass, both signs of `B.n`
each appear more than 10,000 times, and all 256 rotation and polarity pairs
remain solved with compression invariant within `2e-10` relative error.  The
maximum independent normalized continuity error is finite and no greater than
`1e-10`.  `SHK07` follows `SHK06` in priority and in `SMOKE`.

### SHK08 independent tangential electric-field conservation

**What is tested.** `SHK08` checks both independent tangential components of
the ideal-MHD electric field for every solved SHK15 random state, all twelve
named high-precision SHK05 oblique fixtures, 256 proper-rotation pairs, and 256
magnetic-polarity-reversal pairs.  It also requires transformed states to
retain solved compression.

**Why it is tested.** A norm-only or production-residual check can hide a
component permutation, a cross-product sign error, use of heliocentric rather
than shock-frame velocity, or cancellation between tangential components.
Rotation and field reversal make those convention errors observable even when
one special axis-aligned fixture happens to pass.

**How it is tested.** The test reconstructs `u=V-Vsh*n`, expands
`E=-u x B` component-by-component in long-double Cartesian arithmetic, builds
two stable orthonormal tangent axes from each random normal, and projects both
upstream and downstream fields onto each axis separately.  Each component
difference is normalized by its larger physical magnitude, with a `1e-12`
fraction of the full electric magnitude as a near-zero guard.  The calculation
does not call the production cross helper or read `electric_residual`.

**What is expected.** At least 49,000 solved states and exactly twice as many
tangential components are checked; all twelve deterministic fixtures and all
256 transformed pairs participate.  Every transformed state stays solved with
compression invariant within `2e-10` relative error.  Each independent
component residual is finite and no greater than the production threshold
`1e-8`.  `SHK08` follows `SHK07` in priority and in `SMOKE`.

### SHK09 independent momentum-flux conservation

**What is tested.** `SHK09` reconstructs the normal and two tangential
components of the full ideal-MHD momentum flux for every solved SHK15 stress
state and all twelve named SHK05 high-precision oblique fixtures.

**Why it is tested.** Correct mass flux does not guarantee correct momentum
balance.  Pressure, magnetic pressure/tension, shock-frame velocity, or vector
component errors can compensate in one scalar diagnostic while violating a
different tensor component.

**How it is tested.** From serialized primitives the test independently forms
`rho*u_n*u + n*(p+B^2/(2*mu0)) - B_n*B/mu0` in long-double arithmetic and
projects it onto a test-owned orthonormal shock basis.  Each component is
normalized by the larger sum of its dynamic, thermal, and magnetic-stress
magnitudes on either side, avoiding false agreement caused by cancellation.
The stored production `momentum_residual` and production vector helpers are not
used.  A failure reports the component and all three physical term scales with
the complete reproducible fixture record.

**What is expected.** At least 49,000 random solved states and all twelve
deterministic reference states participate.  Normal and both tangential
normalized residuals are finite and individually no greater than `1e-8`.
`SHK09` follows `SHK08` in priority and in `SMOKE`.

### SHK10 independent total-energy-flux conservation

**What is tested.** `SHK10` independently evaluates total ideal-MHD energy
flux for every solved SHK15 state and all twelve named SHK05 fixtures, spanning
the campaign's gamma, beta, Mach-number, obliquity, field, density, and
tangential-flow ranges.  The frozen fixtures must also retain their recorded
sub-`1e-50` high-precision reference residuals.

**Why it is tested.** Energy is the most coupled jump condition: pressure
closure, gamma, all velocity components, magnetic energy, and magnetic work
enter simultaneously.  A missing or duplicated contribution can remain hidden
when compression and simpler invariants are inspected separately.

**How it is tested.** From serialized states, long-double test code separately
forms kinetic transport, enthalpy transport, magnetic-energy advection, and
magnetic-work terms in the shock frame.  Upstream/downstream totals are
normalized by the larger sum of absolute physical contributions, preventing
term cancellation from creating an artificially small denominator.  Neither
the production energy helper nor `JumpResult.energy_residual` is used for
acceptance.  Failures report every term scale and the full fixture record.

**What is expected.** At least 49,000 random states and all twelve independent
fixtures participate, all stored high-precision references remain below
`1e-50`, and every independently normalized energy-flux residual is finite and
no greater than `1e-8`.  `SHK10` follows `SHK09` in priority and in `SMOKE`.

### SHK11 comprehensive shock physical admissibility

**What is tested.** `SHK11` independently classifies all 100,000 SHK15
outcomes and the twelve SHK05 reference solutions.  Solved states are checked
for positive primitives, compression within the gamma-dependent strong-shock
bound, entropy increase, upstream super-fast flow, downstream sub-fast but
super-normal-Alfven flow, and the evolutionary-fast branch identity.  Physical
no-shock and documented weak/singular numerical-limit states have separate
contracts.

**Why it is tested.** Small conservation residuals alone do not identify the
physical fast-shock root.  Intermediate or switch roots can conserve the same
fluxes, and a rejected numerical root must not be silently clipped or exposed
as a converged state to SEP source physics.

**How it is tested.** Characteristic speeds are recomputed with the independent
test-side fast-mode formula; normal flow, normal Alfven speed, entropy proxy,
and compression bound are reconstructed directly from serialized primitives.
No-shock results must preserve the upstream record exactly.  Supported-limit
rejections must retain finite diagnostics while setting `solver_converged`
false.  SHK05 metadata must report exactly one physical root, multiple
converged seeds, and the `EVOLUTIONARY_FAST` branch.

**What is expected.** Every stress case satisfies exactly its applicable
contract, all twelve named references are uniquely evolutionary-fast, and no
generic no-bracket, invalid-state, conservation, or wrong-branch rejection is
present.  The population includes at least 49,000 solved states, exactly
25,000 no-shock states, and at least 25,000 explicit numerical-limit states.
`SHK11` follows `SHK10` in priority and in `SMOKE`.

### SHK17 shock-reference regeneration

**What is tested.** `SHK17` verifies reproducible generation and production
independence of the SHK05 and SHK12 high-precision fixtures.  It checks the
reviewed generator and fixture hashes, byte-for-byte canonical regeneration,
recorded solver/environment metadata, representative higher-precision repeats,
and generator dependencies.

**Why it is tested.** A numerical fixture is not independent evidence if a
production change can silently regenerate it, if its tool environment is
unknown, or if it imports the same reduced equations or branch selector being
validated.  Rounded reference values must also remain stable when reference
precision is increased.

**How it is tested.** `reference/shock_reference_manifest_v1.json` pins the
standard-library Decimal solver version, Python requirement, canonical 80-digit
and audit 100-digit precision, generator hashes, fixture hashes, case counts,
and convergence/branch metadata.  `reference/verify_shock_references.py` runs
both raw-grid generators in an isolated subprocess, compares captured bytes
with the checked-in headers, repeats selected strong, polarity, and weak
families at 100 digits, and requires identical emitted binary64 physics
literals.  AST/token analysis rejects production shock imports or routine
calls.  The verifier never overwrites a checked-in fixture.

**What is expected.** Both regenerated headers are byte-identical to their
canonical SHA-256 values, all selected 100-digit results round to the same
binary64 physics literals as the 80-digit calculations, and the production-
dependency guard passes.  Any generator, fixture, environment, branch, or
conditioning change requires an explicit reviewed manifest/version update.
`SHK17` follows `SHK11` in priority and in `SMOKE`.

Typical direct use is:

```sh
./output/test_swcme --test SHK01
./output/test_swcme --test SHK04
./output/test_swcme --test SHK05
./output/test_swcme --test SHK12
./output/test_swcme --test SHK10
./output/test_swcme --test SHK13
./output/test_swcme --test SHK14
./output/test_swcme --all
```

A shock test failure must not be repaired by increasing the compression floor
or loosening conservation tolerances.  Diagnose the frame transformation,
normal direction, upstream state, root bracket, and downstream reconstruction
first.  A super-fast state for which the nonlinear solver cannot identify an
admissible root is reported with `solver_converged=false` and must not be used
for SEP source physics.


## CON01-CON10: observer-shock magnetic connectivity and cobpoint tracking

`CON01`-`CON10` validate the production 3-D connectivity API implemented by
`swcme3d::Model::observer_connectivity()`.  The connectivity solver does not
maintain a second copy of the shock model: candidate points are tested against
`shape_radius_normal()` and final cobpoints obtain their local plasma/shock
state from `shock_state_direction()`.

The observer Parker line is analytical.  For the same Parker field used by the
production field evaluator,

```text
Delta phi(r) = -Omega (r-r_obs) / V_sw.
```

The observer radial direction is rotated around the configured solar axis by
this amount; colatitude remains constant.  The production parameter
`solar_rotation_rate_rad_s` is shared by field evaluation and connectivity, so
there cannot be a hidden Omega mismatch between the two modules.  The default
value remains the SWCME solar-rotation convention; zero rotation is allowed to
exercise the exact radial reference case.

The intersection residual is

```text
h(r) = r - R_shock[u_Parker(r)].
```

The solver scans from a configurable inner radius to the observer.  It refines
sign-changing roots with bisection, searches local minima of `abs(h)` so tangent
roots that merely touch zero are not missed, and explicitly refines transitions
of the finite-SSE surface-existence flag.  All accepted roots must satisfy the
configured surface-residual tolerance.  Roots are retained in increasing
radius; the default selected cobpoint is the outermost root, corresponding to
the first surface encountered while tracing inward from the observer.

Each `ConnectivityRoot` stores:

- cobpoint radius and Cartesian position;
- signed surface residual;
- analytical Parker arc length from cobpoint to observer; and
- the complete production `LocalShockState`, including `has_shock`, normal,
  normal speed, `theta_Bn`, fast Mach number, compression, and full upstream /
  downstream MHD states.

A geometrical connection and a physical fast shock remain separate concepts.
The connectivity state can therefore be geometrically connected while the
embedded `LocalShockState::has_shock` is false; an SEP source must check the
latter before injection.

The tests are:

- `CON01` — zero-solar-rotation radial limit.  A radial observer line is
  intersected with sphere and SSE geometries and compared with exact analytic
  radius/position/path-length references.
- `CON02` — nonzero-rotation Parker line intersecting a Sun-centered sphere.
  The sphere fixes the root radius exactly while an independent Rodrigues
  rotation verifies the longitude/sign of the production Parker mapping.
- `CON03` — no connection to a finite-width SSE shock.  A field line wholly
  outside the configured cap must return `connected=false` and no fabricated
  flank root.
- `CON04` — exact and near-tangent SSE connection.  The exact half-width ray is
  retained as a valid tangent root, a slightly interior ray connects, and a
  slightly exterior ray remains disconnected.
- `CON05` — multiple intersections and deterministic root selection.  A tightly
  wound Parker line through a strongly non-spherical ellipsoid generates
  multiple physical roots; all are returned in radial order and the outermost
  selection is stable under scan refinement.
- `CON06` — time-continuous cobpoint history.  A ballistic finite SSE front
  evolves from disconnected to connected; the history must show one physical
  onset and continuous/monotone cobpoint motion thereafter without one-step
  classification flicker.
- `CON07` — cobpoint-to-`ShockState` consistency.  Shock quantities stored in a
  connectivity root must be identical to a direct production
  `shock_state_direction()` query at that cobpoint direction.
- `CON08` — analytical Parker path length.  The path length carried by a
  cobpoint is compared with an independent closed-form arc-length reference
  and, away from the pole, must exceed simple radial separation.
- `CON09` — connectivity resolution-limit contract.  Below-, exact-, and
  above-budget requests, an automatic Parker-phase overrun, a maximum-size
  request, a narrow analytic connection window, and the SEP adapter are checked
  to ensure insufficient resolution cannot masquerade as a physical result.
- `CON10` — observer-domain classification.  Null, origin, sub-domain,
  non-finite, exact-boundary, one-ULP-above, negative-axis, and rotated observer
  positions are classified before tracing, including history and SEP-adapter
  paths.

### CON09 detailed validation procedure

What is tested: the test exercises the public
`CONNECTIVITY_SCAN_INTERVAL_BUDGET` boundary and all observable resolution
diagnostics in `ConnectivityState`.  It covers a caller-requested interval
count below the bound, exactly at the bound, one interval above the bound, and
`std::numeric_limits<std::size_t>::max()`.  It separately covers an automatic
Parker-phase sampling requirement that exceeds the bound even though the
caller's radial request is small.  Finally, it checks the
`Interface3D::source_at_observer_cobpoint()` integration path.

Why it is tested: a finite SSE cap can intersect a Parker line over a radial
window narrower than one scan cell.  The former implementation silently
replaced every effective request above 200,000 with 200,000 and then returned
ordinary `Connected` or `Disconnected`.  A missed short window could therefore
be recorded as a scientific conclusion even though the requested numerical
resolution was never run.  The adapter could compound that ambiguity by
mapping every non-connected state to `NoConnection`.

How it is tested: below and at the limit, a zero-rotation Sun-centered sphere
provides an analytic one-root reference and the test runs the actual production
scans, requiring requested and achieved counts to match.  Above the limit, it
requires reject-before-scan diagnostics and no roots.  The phase-derived case
sets the rotation rate so the independently calculated phase request is
`budget + 0.25` steps, whose ceiling is unambiguously `budget + 1`.  The
maximum-size case checks overflow-safe preflight.  For the narrow-window case,
an independent Rodrigues rotation aligns a one-microradian SSE apex with the
Parker line at an off-grid shock radius and verifies the analytic zero surface
residual before submitting an over-budget request.  The same request is then
sent through the SEP adapter.

What is expected: complete below/at-budget scans return the known connected
answer with `requested_scan_intervals == achieved_scan_intervals` and the
published budget.  Every over-budget path returns
`ConnectivityStatus::ResolutionLimit`, `connected == false`, an empty root
list, zero achieved intervals, and the effective requested count plus budget.
The adapter returns `StatusCode::ResolutionLimit`, leaves
`connection_evaluated == false`, and never reports `NoConnection`.  No cap
exhaustion may be reported as unqualified `Connected` or `Disconnected`.

CON09 follows KIN09 in the priority-ordered `SMOKE` profile.  `ROUTINE`, `FULL`,
and `EVENT` include it through `@ALL`.

### CON10 detailed validation procedure

What is tested: CON10 covers the observer-coordinate contract of
`Model::observer_connectivity()`, the stationary-observer history builder, and
`Interface3D::source_at_observer_cobpoint()`.  It checks a null pointer, the
origin, positive- and negative-axis radii one ULP below the lower domain,
individual NaN and positive/negative infinity coordinates, the exact boundary,
one ULP above it, and valid vectors containing negative Cartesian components.
It also exercises equal and reversed source-observer radii.

Why it is tested: an invalid observer must not be scientifically interpreted as
a Parker line that was traced successfully but did not intersect the shock.
That distinction is especially important in connectivity-onset histories and
SEP injection campaigns, where silently counting invalid positions as
disconnected biases onset/loss statistics.  Cartesian signs must not be used as
a proxy for heliocentric radius validity, and the inclusive lower boundary must
not be lost through a separate search-ordering check.

How it is tested: `std::nextafter` constructs the exact one-ULP neighbors of
the shared `MIN_RADIUS_M` constant.  Invalid cases require
`InvalidObserver`, zero requested/achieved scan intervals, and no roots; one
case combines an invalid observer with invalid options to prove observer
validation occurs first.  A zero-rotation spherical fixture gives an analytic
ordinary classification for valid rotated positions.  At the exact lower
boundary, the test evaluates the equal-radius one-point interval twice: first
with a shock elsewhere, then with a stationary spherical shock exactly at the
observer to prove the solver can return both Disconnected and Connected rather
than using a hard-coded boundary result.  Invalid history samples must retain
the same status without tracing.  Adapter calls independently exercise finite
sub-domain, non-finite, reversed-interval, and valid-disconnected cases.

What is expected: null, non-finite, origin, and sub-domain observers return
`ConnectivityStatus::InvalidObserver` before scan or geometry work.  Exact and
above-boundary observers return ordinary `Connected` or `Disconnected` with
finite observer and scan diagnostics; a shock coincident with the exact
boundary produces one zero-path-length root.  Negative coordinate components
remain valid when the vector norm is valid.  The SEP adapter maps finite
sub-domain input to `OutsideModelDomain`, non-finite input to `NonFiniteInput`,
and bad connectivity options to `InvalidConfiguration`, always with
`connection_evaluated=false`.  Only a completed physical disconnection maps to
`NoConnection` with `connection_evaluated=true`.

CON10 follows CON09 in the priority-ordered `SMOKE` profile.  `ROUTINE`, `FULL`,
and `EVENT` include it through `@ALL`.

The default connectivity search tolerance is much tighter than the validation
campaign's `1e-8 AU` root-position target.  Tests deliberately repeat selected
cases with different radial scan densities so correctness cannot depend on a
particular subdivision.  Tolerances must not be relaxed merely to mask missed
roots or unstable tangent classification.

Typical direct use is:

```sh
./output/test_swcme --test CON01
./output/test_swcme --test CON04
./output/test_swcme --test CON06
./output/test_swcme --test CON08
./output/test_swcme --test CON09
./output/test_swcme --test CON10
./output/test_swcme --all
```


## 1D3D01-1D3D03: common-core dimensional-equivalence tests

These tests qualify the common-core refactor rather than a new physical model.
They deliberately configure the 3-D model as a Sun-centered sphere propagating
along +X with the solar rotation axis along +Z and compare it with an equatorial
1-D ray.  In this limit the geometry is exactly equivalent, so any difference is
a software duplication/regression rather than a legitimate dimensional effect.

- `1D3D01` compares the canonical `StepState::common` solar-wind cache and apex
  kinematics field-by-field, then evaluates both public model interfaces at
  1 AU.  Density, radial velocity, `Br`, and `Bphi`/Cartesian azimuthal field
  must agree to roundoff; transverse 3-D velocity and meridional field must
  vanish in the chosen symmetry plane.
- `1D3D02` compares the complete shock state at the spherical apex.  The 1-D
  `StepState::shock_jump` is compared against the 3-D
  `shock_state_direction()` result for shock existence, solver status, radius,
  normal speed, `theta_Bn`, fast-mode speed/Mach number, compression, upstream
  density, full upstream/downstream velocity and magnetic-field vectors,
  downstream density, and downstream pressure.

The tolerances are roundoff-level because both interfaces are now expected to
call the same production core.  A failure must be diagnosed as a wrapper input
mismatch, a reintroduced duplicate equation, or a geometry reduction error; the
tolerance must not be relaxed to hide the difference.

Typical use:

```sh
./output/test_swcme --test 1D3D01
./output/test_swcme --test 1D3D02
./output/test_swcme --all
```

`1D3D03` now validates the shared SOURCE acceleration record.  In the exact
spherical/+X reduction it requires 1-D and 3-D source position, normal, shock
speed, compression, `theta_Bn`, Mach number, upstream density/|B|, DSA
phase-space slope, relative area weight, and deterministic serialized record to
be identical before transport is allowed to diverge.

### Shock-state query-location policy

A local shock state is defined by `(time, shock-surface direction)` and not by
the radius of a background-field sample.  The production sequence is:

1. intersect the selected geometry along the requested direction;
2. evaluate upstream Leblanc/Parker plasma at that surface point;
3. compute the self-similar normal shock speed at that surface;
4. solve the shared ideal-MHD jump;
5. let field/mesh/connectivity consumers use that immutable result.

This ordering is physically important because density and magnetic field vary
with heliocentric radius.  Sampling them at an arbitrary query point would make
Mach number and compression depend on the observer's diagnostic location rather
than on the shock.  `SHK13` and `SHK14` are permanent regression tests for this
contract.

## REG01-REG05: sheath/ejecta region model and mode validation

`REG01`-`REG05` qualify the shared `swcme_regions.hpp` contract used by the 1-D
and 3-D field evaluators.  The exact RH discontinuity remains available in the
shock diagnostic API, while transport-facing FULL_ICME fields use the common C1
shock layer selected by RESOLVED_COMPRESSION.  SHOCK_ONLY/SOURCE keeps the
transport background analytical on both sides of the mathematical source.

- `REG01` — **SHOCK_ONLY upstream-field identity**.  Samples points ahead of
  and geometrically behind the expanding front in both 1-D and 3-D.  Density,
  velocity, and Parker magnetic field must be identical to the analytical
  upstream state.  The test deliberately changes FULL_ICME sheath/ejecta
  parameters to extreme valid values and proves they have no effect in
  SHOCK_ONLY mode.
- `REG02` — **FULL_ICME resolved-shock inner RH boundary**.  Evaluates the
  inner endpoint `R_sh-w_sh/2` of the RESOLVED_COMPRESSION C1 shock layer and
  verifies that density, velocity, and magnetic field equal the complete
  production Rankine-Hugoniot downstream state.  This preserves the exact RH
  boundary while allowing the transport field to resolve the jump numerically.
- `REG03` — **magnetic-ejecta density and velocity factors**.  Evaluates the
  middle of the ejecta, away from LE/TE blends, for factors below, equal to,
  and above unity.  `f_ME=0.5` must give `0.5*n_up` and
  `V_ME_factor=0.8` must give `0.8*V_sw`.  Negative factors must fail
  centralized configuration validation rather than be clipped.
- `REG04` — **self-similar local layer nesting**.  Samples a finite SSE cap from
  apex to near-flank and verifies
  `R_LE/R_sh = 1-f_sheath` and
  `R_TE/R_sh = 1-f_sheath-f_ejecta` at every direction.  Layer ordering must
  never invert.  A public configuration whose thickness fractions sum to one
  or more is rejected.
- `REG05` — **continuity/smoothness at artificial region transitions**.  Checks
  the common symmetric LE/TE smoothstep weights and samples both public model
  interfaces around every transition endpoint.  Equivalent spherical +X
  1-D/3-D fixtures must agree to roundoff.  One-sided finite-difference
  derivatives converge across the C1 artificial interfaces.  The resolved
  shock itself is covered separately by `ACC03`/`ACC04` because its smoothing
  is controlled by the acceleration representation rather than by LE/TE region
  phenomenology.

Run them directly with:

```sh
./output/test_swcme --test REG01
./output/test_swcme --test REG02
./output/test_swcme --test REG03
./output/test_swcme --test REG04
./output/test_swcme --test REG05
./output/test_swcme --all
```

The region tests must not be made green by reintroducing a compression floor,
clipping sub-unity ejecta factors, or sorting invalid boundary radii at runtime.
Such behavior defeats the physical/configuration contracts the tests are meant
to protect.


## ACC01-ACC05: single shock-acceleration representation

These tests qualify `swcme_acceleration.hpp` and the shock-smoothing portion of
`swcme_regions.hpp`.  Their purpose is to prevent one SEP population from seeing
a pre-imposed DSA source and a second resolved compression accelerator at the
same shock.

- `ACC01` — **SOURCE explicit source / SHOCK_ONLY flow**.  Requires a physical
  fast shock to produce `source_enabled=true` and
  `resolved_compression_enabled=false`, checks the phase-space DSA slope
  `q=3r/(r-1)`, and samples both sides of the mathematical front to prove the
  transport background remains the exact Parker/Leblanc solar wind.
- `ACC02` — **mutual-exclusion validation**.  Rejects `SOURCE+FULL_ICME`,
  `RESOLVED_COMPRESSION+SHOCK_ONLY`, zero shock width in resolved mode, and a
  negative relative source weight.  These are configuration errors rather than
  runtime conventions.
- `ACC03` — **resolved-compression shock profile**.  Verifies a finite total
  width, exact upstream outer endpoint, exact RH downstream inner endpoint,
  an intermediate midpoint, and zero-slope/C1 matching at both ends of the
  common numerical shock layer.
- `ACC04` — **1-D/3-D resolved-profile identity**.  Samples five normalized
  positions across an equivalent spherical/+X shock and requires density,
  radial velocity, Br/Bx, and Bphi/By to agree to roundoff.  This prevents
  dimension-dependent shock smoothing from masquerading as a transport effect.
- `ACC05` — **resolved mode disables prescribed DSA source**.  Requires
  `source_enabled=false`, `resolved_compression_enabled=true`, no active DSA
  slope, zero source weight, and an explicit `NA` slope in deterministic audit
  serialization.

`1D3D03` complements these tests by comparing the complete SOURCE acceleration
record from equivalent 1-D and 3-D configurations and requiring byte-identical
serialization.

Run the acceleration gates directly with:

```sh
./output/test_swcme --test ACC01
./output/test_swcme --test ACC02
./output/test_swcme --test ACC03
./output/test_swcme --test ACC04
./output/test_swcme --test ACC05
./output/test_swcme --test 1D3D03
```

A passing test must not be obtained by merely masking `div(V)` after constructing
an RH velocity jump in SOURCE mode.  SOURCE removes the resolved shock from the
transport background by using SHOCK_ONLY; RESOLVED_COMPRESSION owns the single
finite-width compression profile.

## SEP01-SEP06: SWCME-to-SEP/AMPS integration contract

These tests qualify `swcme_sep_source.hpp` and `swcme_sep_interface.hpp`.  The
integration layer is a consumer of the production SWCME model, not a second
shock/connectivity implementation.

- `SEP01` — **background adapter identity**.  A 3-D AMPS-facing single-point
  query must reproduce the direct checked production density, velocity,
  magnetic field, `|B|`, and `div(V)` result in SI.
- `SEP02` — **spectrum/unit convention**.  Checks differential-intensity unit
  conversion round trips, relativistic proton rigidity, DSA slope-to-intensity
  indices, unity at `E_ref`, and physical `J(E_ref)` when explicit reference
  normalization is selected.
- `SEP03` — **1-D/3-D complete source-record identity**.  Equivalent spherical
  +X configurations must emit byte-identical deterministic `SEPSourceState`
  records and matching adapter background values.  This is the final interface
  guard beyond lower-level `1D3D03`.
- `SEP04` — **finite-SSE source-surface weighting**.  Builds source records from
  the corrected triangular mesh, requires physical patch areas, normalizes
  active area fractions to unity, and verifies that total relative patch
  weight equals the configured uniform source weight.
- `SEP05` — **cobpoint source reuse**.  Obtains an observer source through the
  production Parker connectivity solver and verifies that source position and
  shock quantities correspond to the selected production cobpoint rather than
  a separately solved surface.
- `SEP06` — **resolved-compression/source exclusion and manifest stability**.
  `RESOLVED_COMPRESSION` must expose no prescribed source spectrum,
  `relative_intensity_shape()` must return `SOURCE_INACTIVE`, and the resolved
  SWCME+SEP manifest must be deterministic.

`SEPSourceState` makes source and connectivity state explicit through
`active`, `connection_evaluated`, and `connected`.  Its geometry and shock
fields use SI units.  Default `RELATIVE_ONLY` normalization remains
dimensionless; a dimensional spectrum exists only when
`REFERENCE_DIFFERENTIAL_INTENSITY` supplies a positive SI `J(E_ref)`.

Run the integration gates directly with:

```sh
./output/test_swcme --test SEP01
./output/test_swcme --test SEP02
./output/test_swcme --test SEP03
./output/test_swcme --test SEP04
./output/test_swcme --test SEP05
./output/test_swcme --test SEP06
```

## ERR01-ERR05: explicit numerical-status propagation

These tests enforce the rule that the physics layer may not convert a failed
numerical state into a plausible physical value.  They qualify the shared
`swcme_status.hpp` status contract, the checked 1-D/3-D evaluators, the explicit
RH solver outcome, and output-data validation.

- `ERR01` — **1-D outside-domain query**.  Requests a radius below the supported
  `1.05 R_sun` Parker/Leblanc boundary and requires
  `OUTSIDE_MODEL_DOMAIN` with `sample_index=0`.  Sentinel outputs at the failed
  sample must remain unchanged.  The source-compatible void evaluator must
  throw rather than clipping the radius to the domain boundary.
- `ERR02` — **3-D non-finite Cartesian input**.  Inserts a NaN into the second
  point of a two-sample batch and requires `NONFINITE_INPUT` with
  `sample_index=1`.  The failed sample must not be rewritten as zero/ambient
  plasma, and the legacy wrapper must throw.
- `ERR03` — **Rankine-Hugoniot outcome classification**.  A degenerate normal
  and a non-finite upstream magnetic component must be `INVALID_INPUT`; a
  well-formed sub-fast front must be `NO_SHOCK`; a well-conditioned fast shock
  must be `SOLVED`.  The test prevents a bad primitive state from masquerading
  as an unmagnetized/no-shock solution.
- `ERR04` — **no arbitrary +X normalization fallback**.  A zero shock direction
  must return `DEGENERATE_VECTOR`, leave `surface_exists=false`, and make the
  old bool shock wrapper throw rather than constructing a +X surface.
- `ERR05` — **writer rejects corrupt physics**.  A shock mesh containing NaN is
  passed to the checked surface writer.  The writer must return
  `NONFINITE_RESULT` before creating the requested file; substituting zero for
  the bad coordinate is forbidden.

Run the status gates directly with:

```sh
./output/test_swcme --test ERR01
./output/test_swcme --test ERR02
./output/test_swcme --test ERR03
./output/test_swcme --test ERR04
./output/test_swcme --test ERR05
./output/test_swcme --all
```

`NO_SURFACE` and RH `NO_SHOCK` are deliberately not treated as numerical
failures: they represent valid physical/geometrical outcomes.  The SEP layer
similarly treats `NO_CONNECTION` and `SOURCE_INACTIVE` as explicit expected
absence states.  Conversely, a solver, normalization, domain, or
non-finite-value failure must remain visible to the caller and must never be
made green by restoring `finite_or`, radius clipping, or ambient/zero
substitution.

## DEN04: pressure and sound-speed closure

**What is tested.** The test covers both `PROTON_ONLY` and `MULTI_SPECIES`
prepared thermodynamics at gamma values 1.2, 1.4, 1.5, and 5/3. Fixtures include
equal and unequal species temperatures, zero and finite alpha abundance, and a
cold but positive pressure. It also checks the upstream primitive passed by
the public 1-D and 3-D shock APIs and rejects unknown closure values, negative
abundance, non-finite electron temperature, and non-positive alpha temperature.

**Why it is tested.** Pressure and mass density jointly determine sound speed,
fast-mode Mach number, and the MHD jump. A partially applied composition model
could make adapters, dimensional models, and the shock solver describe
different plasma while all individual values remain finite.

**How it is tested.** Test-owned long-double equations use independently
pinned Boltzmann, proton-mass, and alpha-mass literals. Charge neutrality is
solved directly, partial pressures are summed independently, and sound speed is
reconstructed as `sqrt(gamma*p/rho)`. Production values come only from public
prepared state and shock paths.

**Expected result.** Every density, mass, pressure, and sound-speed residual is
below `1e-13`; both model APIs use the selected closure exactly; every invalid
configuration is rejected before preparation. `DEN04` follows `SHK17` in the
priority registry and `SMOKE` profile.

## DEN03: composition and mass density

**What is tested.** `DEN03` verifies electron, proton, and alpha number
densities, charge neutrality, alpha abundance, and total ion mass density for
`f_alpha=0`, 0.04, and 0.10 over `n_e=1e2` through `1e14 m^-3`.

**Why it is tested.** Leblanc supplies electron density, whereas Alfven and
fast-mode speeds require mass density. Treating `n_e` as proton density when
alpha abundance is nonzero biases shock strength without changing displayed
electron density.

**How it is tested.** Test-owned long-double algebra uses
`f_alpha=n_alpha/n_proton`, solves `n_e=n_proton+2*n_alpha`, and combines
independently pinned proton and alpha masses. It checks all fields, the charge
balance, abundance monotonicity, and the exact zero-alpha limit.

**Expected result.** All relative errors are below `1e-13`, charge neutrality
closes, mass density grows with alpha abundance, and the no-alpha result is
exactly `rho=m_proton*n_e`. `DEN03` precedes `DEN04` in the priority registry
and `SMOKE` profile.

## DEN02: Leblanc density asymptotic behavior

**What is tested.** `DEN02` checks every individual Leblanc power-law term and
their sum at 41 logarithmically spaced radii spanning 0.01--100 AU, plus the
required 0.5, 1, 2, and 5 AU checkpoints. It checks both dimensional public
evaluators and the far-field `r^2 n(r)` limit.

**Why it is tested.** Exact normalization at 1 AU cannot detect a wrong radial
power or a compensating coefficient error. The asymptotic test isolates those
errors and demonstrates where the profile has actually entered its `r^-2`
regime.

**How it is tested.** A test-owned long-double oracle pins the published A, B,
and C coefficients and the adopted AU/solar-radius conversions, independently
derives the normalization, and evaluates all three terms. The observed excess
of `r^2 n` above its asymptote is compared with `(n4+n6)/n2` and required to
decrease monotonically.

**Expected result.** Each term and total-density residual is below `1e-12`.
The scaled density remains above and approaches its asymptote, with its excess
fully explained by the higher-order terms. `DEN02` follows `DEN04` in priority
and in `SMOKE`.

## DEN05: default closure compatibility

**What is tested.** A default configuration is compared with an explicitly
selected `PROTON_ONLY` configuration through 1-D fields, 3-D fields, both shock
paths, AMPS background pressure/focusing, serialized SEP source records, and
resolved configuration provenance.

**Why it is tested.** Introducing an optional composition closure must not
silently change established campaigns or make old default construction select
new physics.

**How it is tested.** Independent model and adapter instances are prepared at
the same time. Every floating-point output is compared bit-for-bit (owner IDs
are correctly excluded), and source CSV records are compared byte-for-byte.
The test separately opts into `MULTI_SPECIES` and checks its manifest name and
distinct configuration digest.

**Expected result.** Default and explicit proton-only results are exact; the
multi-species mode appears only after explicit selection and has distinct
provenance. `DEN05` follows `DEN02` in priority and in `SMOKE`.

## PAR04: Parker field solenoidality

**What is tested.** Cartesian `div(B)` is measured at 1,024 deterministic
random positions with random solar axes, wind speeds, rotation rates, Parker
source radii, normalization latitudes, magnetic strengths, and both global
polarities. Each point uses four successively halved relative stencil widths.

**Why it is tested.** Component and magnitude benchmarks can pass while the
assembled 3-D vector violates the Maxwell solenoidal constraint. Axis rotation
and a finite source radius expose basis and radial-pitch defects hidden by an
equatorial default fixture.

**How it is tested.** The test computes centered Cartesian derivatives by
calling the production field at six independently perturbed points per stencil.
It normalizes `|div B|` by `|B|/r`, selects the truncation/roundoff plateau from
four widths, checks second-order improvement, and repeats the finest estimate
after exact global polarity reversal.

**Expected result.** The median plateau residual is below `1e-8`, the maximum
below `1e-6`, at least 90% of points demonstrate second-order improvement, and
polarity changes residuals by less than `1e-12`. `PAR04` follows `DEN05` in
priority and in `SMOKE`.

## PAR05: Parker field-line tangency

**What is tested.** The analytical observer-anchored Parker line, Cartesian
background field, accumulated connectivity phase, and finite-source-radius
pitch are compared in 1,024 random observer/radius/configuration cases. Both
global magnetic polarities are covered.

**Why it is tested.** A field and connectivity map can each look plausible but
still use different signs, rotation axes, or source-radius offsets. That defect
would place a cobpoint off the actual field line and corrupt field-aligned SEP
transport distances.

**How it is tested.** Test-owned Rodrigues rotation and the independently
integrated `dphi/dr=-Omega(1-rb/r)/V` construct the reference curve. Centered
radial differentiation produces its tangent; the production magnetic vector is
queried at that point. The test also projects `Br` and `Bphi` onto independently
constructed bases and checks their ratio against the analytic pitch.

**Expected result.** Connectivity positions agree within `2e-13` relative,
pitch differs by less than `2e-13`, and the tangent-field angle is below
`1e-8` rad—parallel for positive polarity and antiparallel for negative.
`PAR05` follows `PAR04` in priority and in `SMOKE`.

## PAR06: Parker path length and focusing

**What is tested.** Closed-form Parker arc length, magnetic focusing length,
their SI units and signs, finite source radius, and the exact zero-rotation
limit are exercised in 160 randomized fixtures.

**Why it is tested.** Focused SEP transport uses distance along the field and
`d ln|B|/ds`; a vector field can pass pointwise checks while these derived
transport quantities retain a sign, offset, or differentiation defect.

**How it is tested.** A test-owned adaptive Simpson integrator evaluates
`sqrt(1+k^2(r-rb)^2)`. A five-point, fourth-order stencil independently
differentiates `ln|B|` along radius, converts it to the outward path coordinate,
and is checked at two decreasing step sizes before comparison.

**Expected result.** Relative path error is below `1e-9`, focusing error below
`1e-7`, the zero-rotation length equals the radial interval exactly, and
baseline outward focusing is finite and positive. `PAR06` follows `PAR05` in
priority and in `SMOKE`.

### PAR02 expanded random-vector campaign

**What is tested.** The original transparent latitude/axis fixtures are joined
by 1,024 fixed-seed vectors spanning radius, longitude, latitude, hemisphere,
solar axis, field polarity, wind speed, rotation rate, source radius, and 1-AU
normalization latitude.

**Why it is tested.** Equatorial examples do not fully exercise basis
orientation, hemispheric behavior, source-radius pitch, or rotational
covariance.

**How it is tested.** The validation constructs `e_r` and `e_phi` independently,
normalizes signed `Br` at 1 AU, computes the source-aware `Bphi`, and transforms
to Cartesian. Paired models reverse polarity, while a Rodrigues rotation is
applied independently to the axis, sample point, and reference vector.

**Expected result.** Vector and rotation-covariance errors are below `1e-12`;
polarity reverses the complete vector to roundoff without changing magnitude.

### PAR03 expanded polar-limit campaign

**What is tested.** Both Parker poles are approached for four arbitrary axes,
both radial polarities, a finite source radius, eight longitudes, and six
logarithmic angular offsets from `1e-2` to `1e-12` rad; exact poles are included.

**Why it is tested.** `Omega_hat cross e_r` vanishes at a pole. Normalizing it
without a regular limiting branch can amplify roundoff into a finite,
longitude-dependent transverse field.

**How it is tested.** An independent tangent basis constructs every approach.
The returned field is projected onto `e_r` and `e_phi`; `Bphi` is compared with
the source-aware analytic `sin(theta)` scaling, and equal-offset magnitudes are
compared across longitudes.

**Expected result.** All vectors remain finite, exact poles are radial,
transverse magnitude vanishes as `O(sin(theta))`, and no longitude-dependent
limiting residue appears.

### PAR01 expanded equatorial component campaign

**What is tested.** Equatorial radial and azimuthal Parker components, their
Cartesian orientation, `r^-2` radial scaling, finite-source-radius pitch, and
both global magnetic polarities are exercised in 96 matrix cases. The original
+X point is retained, and the canonical-axis matrix explicitly includes +X and
+Y along with an arbitrary proper rotation axis.

**Why it is tested.** A vector can have the right magnitude while carrying an
incorrect azimuthal sign, polarity convention, radial scaling, or source-radius
offset. Those defects directly reverse field-line transport or move modeled
magnetic footpoints.

**How it is tested.** Each point is independently decomposed in the equatorial
`(e_r,e_phi)` basis. A test-owned analytical Parker construction supplies the
Cartesian reference; separate projections check the sign of `Br` and the
relation `Bphi/Br=-Omega(r-rb)/Vsw` for two axes, two directions, two
polarities, three radii, two wind speeds, and two source radii.

**Expected result.** Every Cartesian vector and projected pitch agrees within
`1e-12`, `Br` has the configured polarity, and the zero-source-radius cases
remain compatible with the historical winding law.

## CROSS01: one-dimensional and three-dimensional physics consistency

**What is tested.** Equivalent radial 1-D and spherical 3-D configurations are
compared for density, thermodynamic pressure, velocity, magnetic components and
magnitude, divergence, Parker path and focusing lengths, shock radius and
normal speed, compression, `thetaBn`, Mach number, and acceleration-source
state at both apex and flank directions.

**Why it is tested.** Sharing internal helpers does not prove that two public
wrappers transfer configuration, interpret geometry, and expose statuses in
the same way. Cross-model drift can otherwise survive unit-level tests.

**How it is tested.** Four configurations cover two closure policies, two
polarities, zero/finite source radius, and strong/weak/no-shock behavior. Three
times and three radii yield 36 field comparisons, while two equatorial sphere
directions yield 24 shock/source comparisons. Dimension-aware relative errors
are used for values ranging from Tesla to meters per second.

**Expected result.** All status flags are identical and every finite common
quantity agrees within `3e-12` (most field quantities within `3e-13`).

## GEO01: geometry rotation and derivative convergence

**What is tested.** Sphere, finite-width SSE, and triaxial ellipsoid positions,
normals, radii, normal speeds, compression, and local area elements are tested
under an arbitrary proper rotation. Four centered time steps test normal-speed
derivative convergence.

**Why it is tested.** Hidden dependence on global Cartesian axes changes a
nonspherical CME when coordinates are rotated, while a pointwise geometry test
can miss an incorrect time derivative used in shock fluxes.

**How it is tested.** Test-owned Rodrigues rotation transforms every physical
vector. Surface results are transformed back or compared as scalar invariants;
local areas come from independent three-point cross products. Centered temporal
differences at 16, 8, 4, and 2 seconds reconstruct surface normal speed.

**Expected result.** Position and normal errors remain below `3e-13`, local
area is invariant within `2e-10`, scalar diagnostics agree within `3e-12`, and
each shape demonstrates at least two convergent derivative refinements with a
final relative error below `1e-7`.

## MESH01: mesh and integrated source convergence

**What is tested.** Triangle validity, analytical-normal orientation,
mean-ratio quality, total surface area, active source-weight normalization, and
an integrated production source measure are checked for sphere, SSE, and
ellipsoid meshes.

**Why it is tested.** A mesh may look plausible while carrying inverted or
poorly conditioned cells, and pointwise source records do not establish that a
surface-integrated injection converges with resolution.

**How it is tested.** Four nested meshes (8x16 through 64x128) are built for
each shape. The test computes triangle quality independently and integrates
`patch_area*Vsh,n*(compression-1)` from actual active SEP source records. The
first three levels are compared with the separately rebuilt finest mesh.

**Expected result.** All cells are finite, positive, and outward; median quality
exceeds 0.05; source weights sum to the configured value; area and integrated
source errors decrease monotonically; final changes are below 0.5% and 1%.

## REG01: region boundary and smoothing convergence

**What is tested.** Classification and weights on both sides of the shock,
leading edge, and trailing edge; zero and finite smoothing; local requested
versus resolved widths; and spatial and temporal derivative convergence.

**Why it is tested.** Boundary inequalities, hidden width clipping, or stale
apex-scaled widths can introduce discontinuities and nonconvergent transport
sources even when representative interior points are correct.

**How it is tested.** Eight logarithmic offsets bracket every boundary. An
independent cubic checks weights and normalization. Centered spatial and time
steps are halved four times, and 64 deterministic SSE position/time samples
verify exact local self-similar widths.

**Expected result.** Zero widths produce explicit region changes without blend;
finite weights match the independent smoothstep, sum to one, and converge in
space. Resolved widths are bit-exact products of the accepted fractions and
local radii, and all boundary velocities converge to their DBM references.

## CON11: connectivity random stress and transition convergence

**What is tested.** Valid-domain classification, complete scan diagnostics,
all-root ordering and residuals, deterministic outer-root selection, analytic
sphere outcomes, onset/loss timing, forward/reverse history equivalence, grid
refinement, and explicit resolution-cap exhaustion.

**Why it is tested.** Sparse hand-picked cobpoints cannot expose rare tangent,
multiple-root, orientation, or state-transition failures, while silently
truncated searches can mislabel an unresolved case as disconnected.

**How it is tested.** A fixed xorshift stream drives 10,000 cases over 100
configurations and 100 observers. Roughly one third use an independent sphere
interval oracle. A moving SSE is sampled on 30- and 15-minute forward grids and
the fine grid in reverse; one request deliberately exceeds the published cap.

**Expected result.** Every ordinary case is Connected or Disconnected with a
complete scan; all roots are finite, ordered, and residual-qualified; analytic
sphere outcomes agree; onset and loss converge within one fine interval;
reverse history is identical; cap exhaustion returns `ResolutionLimit`.

## SAN01: memory and undefined-behavior validation

**What is tested.** The complete deterministic registry, malformed-input paths,
representative and high-count random campaigns, output handling, integration
adapters, and all three demonstration executables are exercised under ASan and
UBSan.

**Why it is tested.** Correct numerical results do not rule out out-of-bounds
access, use-after-free, invalid lifetime behavior, signed overflow, invalid
shifts, alignment faults, or other undefined behavior.

**How it is tested.** `make san01-sanitize` compiles fresh sources with
`-fsanitize=address,undefined -fno-sanitize-recover=all
-fno-omit-frame-pointer`, executes `--all`, and runs each demo in an isolated
directory. A child marker prevents only recursive SAN01 orchestration.

**Expected result.** Compilation, 127 registered validations, and all demos
finish with no sanitizer-origin diagnostic. Leak enumeration is explicitly
disabled by default in ptrace/container environments; set
`SAN01_ASAN_OPTIONS=detect_leaks=1:halt_on_error=1` on a compatible untraced
host to add LeakSanitizer.

## Priorities 50-60: execution, evidence, and observational campaign gates

The registered `V1`-`V5` cases use clearly labelled `SYNTHETIC_REGRESSION`
fixtures to test calculations and acceptance logic. V6 is labelled
`COUPLING_SMOKE_TEST` because it exercises the serialized AMPS boundary and a
validation-owned transport consumer without claiming external AMPS science
skill. Neither fixture class substitutes for an observational data package.
`EVT01` enforces that distinction: a `RELEASE_VALIDATION` campaign cannot pass
without traceable observational provenance and complete evidence for every
V1-V6 layer.

### THR01: thread and scheduler reproducibility

**What is tested.** Prepared-state evaluation under 1, 2, 4, and 8 workers,
static and atomic dynamic work distribution, four chunk sizes, reversed
within-chunk order, repeated interleavings, per-job status/results, work counts,
canonical record hashes, reductions, and unchanged prepared-state seals.

**Why it is tested.** A race, duplicate/missing job, order-dependent reduction,
or mutable cache may appear only under one scheduler even if serial results pass.

**How it is tested.** Ninety-nine heterogeneous jobs are first evaluated by a
serial oracle, then by 48 schedule combinations. Workers write disjoint slots;
the post-join comparison and reduction always use canonical job order.

**Expected result.** Every slot and status is exact, every job is visited once,
hashes match, reductions are bitwise equal or within roundoff, and state seals
remain unchanged. The suite is also compiled and run by SAN01 for race-adjacent
lifetime diagnostics. `make thr01-tsan` builds and runs the same matrix with
ThreadSanitizer on a TSan-capable CI executor.

### COV01: line and branch coverage closure

**What is tested.** Source-only statement and branch coverage, including
invalid input, singular limits, boundary logic, and writer cleanup paths.

**Why it is tested.** Passing numerical fixtures can leave rejection and
cleanup code unexecuted, making regressions in safety behavior hard to detect.

**How it is tested.** `make cov01-coverage` performs a fresh GCC coverage build,
runs a reviewable nonredundant registry set, merges JSON gcov records by source
and line/branch identity, and writes `output/cov01/coverage.json`. Tests,
generated references, and physics exclusions are not counted as production.

**Expected result.** Production coverage exceeds the frozen 72% line and 40%
branch gates. The current verified result is 89.60% (3887/4338) lines and
57.83% (2202/3808) branches; future modified production paths must remain
covered rather than being excluded.

### PERF01: performance and scaling guardrails

**What is tested.** Preparation, scalar/batch evaluation, 32 shock solves,
mesh plus output preflight, integrated source construction, connectivity, and
field batches at 1/2/4/8 threads; median, nearest-rank p95, environment, and an
observable checksum are reported.

**Why it is tested.** Accidental repeated preparation, per-point allocation,
quadratic loops, or lost parallel sharing can preserve answers while making a
campaign impractical.

**How it is tested.** The registered parent builds a fresh `-O3 -DNDEBUG`
standalone child. Each workload is warmed and sampled 7 or 9 times; timing is
explicitly excluded from sanitizer/coverage binaries.

**Expected result.** P95 remains below 100 ms for preparation/batch, 3 s for
shock/mesh, and 5 s for source/connectivity. Each parallel median is less than
2.5 times the serial median and the checksum is finite. These portable CI
ceilings catch severe regressions; controlled release hardware should retain a
separate frozen baseline for tighter trend monitoring.

### REP01: fixture and campaign reproducibility

**What is tested.** Path-ordered input inventory, source/config/event hashes,
seed, selected tests, compiler contract, and the canonical campaign fingerprint.

**Why it is tested.** A result without input/tool identity cannot be recreated
or distinguished from a silently changed fixture.

**How it is tested.** Manifest schema v2 hashes all reference, profile, example,
and Python manager inputs. REP01 builds the same manifest twice and requires
identical fingerprints despite timestamp/host metadata, then changes the seed
and requires a different fingerprint.

**Expected result.** Identical defining inputs reproduce the fingerprint and
inventory byte-for-byte; any defining seed/source/config/fixture/tool-contract
change changes it. Output reference exports retain their own SHA-256.

### EVT01: campaign schema and completeness

**What is tested.** Schema versioning, known keys, provenance, layer presence,
evidence, metrics, uncertainty, convergence, valid NA rationale, and the exact
`PASS`/`FAIL`/`INCOMPLETE`/`NOT_APPLICABLE`/`ERROR` states.

**Why it is tested.** Execution success is not scientific completeness; the old
empty EVENT arrays could incorrectly aggregate to PASS.

**How it is tested.** The Python contract receives complete and deliberately
damaged packages: wrong versions, unknown keys, missing V layers, empty metrics,
synthetic-only release provenance, and a threshold failure. Analysis commands
are not launched when schema state is ERROR or INCOMPLETE.

**Expected result.** Malformed packages are ERROR, missing science is
INCOMPLETE, complete threshold violations are FAIL, justified NA is accepted
only for nonrelease self-tests, and only a complete package can PASS.

### V1: background Parker and Leblanc validation

**What is tested.** The observational metric pipeline for density, magnetic
field magnitude, and Parker angle from 0.30 to 1.40 AU.

**Why it is tested.** Event conclusions require the upstream background to be
credible before CME/shock differences are interpreted.

**How it is tested.** Production 3-D Parker/Leblanc queries are compared with a
fixed synthetic reference perturbation. Ratios and median angular error are
computed independently, and a deliberately bad angle set proves rejection.

**Expected result.** Median angle error is at most 20 degrees and normalized
density/field trends remain within factor two. Real V1 evidence must replace
the synthetic records with quiet-interval observations and traceable hashes.

### V2: CME and shock apex kinematics

**What is tested.** DATA_DRIVEN height knots, DBM propagation, and held-arrival
absolute and relative error metrics.

**Why it is tested.** A correct local shock solver cannot compensate for an
incorrect apex trajectory or biased arrival time.

**How it is tested.** Production PCHIP states are evaluated at every supplied
synthetic knot. A production DBM trajectory is independently bisected at 0.96
AU and compared with a held synthetic arrival offset.

**Expected result.** Knot errors are below `1e-12`; arrival error is at most six
hours and 15%. A real campaign also records the useful 12-hour/25% tier and
uses independently sourced height-time and held arrival observations.

### V3: in-situ shock jump validation

**What is tested.** Compression, downstream speed and field, shock obliquity,
and all five normalized Rankine-Hugoniot residuals.

**Why it is tested.** A plausible arrival does not establish a conservative or
physically admissible local shock jump.

**How it is tested.** The production ideal-MHD solver evaluates a fixed oblique
state; independent metric code compares its result with controlled synthetic
measurement offsets and aggregates the maximum conservation residual.

**Expected result.** Residual is at most `1e-8`; compression, speed, field, and
theta-Bn errors are at most 25%, 15%, 30%, and 15 degrees. Observational PASS
requires actual upstream/downstream windows and uncertainty records.

### V4: shock geometry and encounter validation

**What is tested.** Multi-observer hit/miss classification, apex-to-flank
ordering, arrival timing, and the uncertainty decision.

**Why it is tested.** An apex trajectory may be accurate while width or flank
geometry predicts the wrong spacecraft encounter.

**How it is tested.** Four synthetic observer longitudes bracket a 45-degree
front; a transparent cosine-flank proxy produces ordered arrivals compared
with independent synthetic times and an eight-hour uncertainty.

**Expected result.** All hit/miss decisions and ordering are correct and every
arrival is inside uncertainty. Real evidence must provide ephemerides, local
normal estimates, and documented timing uncertainty.

### V5: magnetic connectivity benchmark

**What is tested.** Qualitative connection histories, transition counting,
four-hour sampling sensitivity, required observers, and multi-event plumbing.

**Why it is tested.** Connectivity is a time-history classification; one
cobpoint or event-specific tuning cannot validate onset and loss behavior.

**How it is tested.** The required 2010-09-09 STA-stays-connected and
STB-connects-then-loses constraints from Tao et al. (2025, ApJ 995:77,
doi:10.3847/1538-4357/ae17bd) seed the qualitative comparator. SOHO and two
additional histories are synthetic, and shifted sampling checks transition
stability without event-specific production tuning.

**Expected result.** History classes and transitions match and survive the
sampling perturbation. The bundled test remains scientifically INCOMPLETE;
release PASS requires traceable STA/STB/SOHO observations for the mandatory
event plus at least two independent observed events and sensitivity evidence.

### V6: SEP-facing AMPS integration

**What is tested.** The complete, serialized SWCME SEP source record is checked
at both standalone-to-adapter boundaries. A separately parsed source history
then drives fixed-source and SWCME-dependent 1-D/3-D transport controls, with
perpendicular diffusion disabled and enabled as distinct attribution cases.

**Why it is tested.** The application depends on a correct SWCME-to-AMPS
handoff. A unit shift, reordered CSV column, hidden adapter normalization, or
uncontrolled random/numerical setting can otherwise be misdiagnosed as either
SWCME physics or AMPS transport skill.

**How it is tested.** `V6` constructs source histories directly from the 1-D
and 3-D production models and independently through their public SEP adapters.
It requires complete record bytes and FNV-1a regression fingerprints to match
for each producer/adapter pair. A validation-owned parser accepts exactly the
37-column public schema and does not share `SEPSourceState`. Parsed records must
remain active, super-fast, monotonically propagating, and carry physical
compression and DSA slope values. A deterministic consumer then holds the
source normalization fixed, activates SWCME source dependence, compares 1-D
with 3-D while perpendicular diffusion is off, and repeats 3-D with it on using
the same `0x56365f414d50535f` seed, 600-second step, and 288-step horizon.

**Expected result.** Both standalone/adapter histories are byte- and
hash-identical; the independently parsed source history is physically
plausible; fixed and source-dependent 1-D/3-D no-perpendicular profiles are
identical; source activation and perpendicular diffusion each produce a
finite, nonnegative, attributable change. The test logs
`coupled_flux_skill=NOT_EVALUATED`: the consumer is a coupling smoke harness,
not the external AMPS solver. A release V6 PASS must additionally archive real
AMPS build/provenance, transport outputs, convergence and uncertainty evidence;
it cannot repair a failed SWCME gate.

Run the focused gate with:

```sh
./output/test_swcme --test V6
```
