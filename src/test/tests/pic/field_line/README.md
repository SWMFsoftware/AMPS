# Field-line geometry verification

The production geometry kernels are in `src/pic/pic_field_line_geometry.h`.
`CORE-PIC-FIELD-LINE-001` runs eleven deterministic cases, including motion of
both endpoints, repeated inner insertion, outer chord/sphere coverage,
failure atomicity, unequal segment exposure and curved magnetic tracing with
either magnetic polarity. Run from the repository root:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-FIELD-LINE-001 --output-dir /tmp/amps-field-line-core
```

The native fixture uses actual `cFieldLine`, particle buffers, segment lists,
sampling and output, with manufactured SI data and fixed dt = 0.75 s. It checks
empty inner segments, retiring particle weights, preserved particle velocity,
shifted material coordinates, ownership-table replacement, shock intersection,
fluence and equal density on segments with unequal ages. Output must record
completed counts 1,1,2,2 and the unchanged timestep. Expected density and fluence
are computed from prescribed weights, volume, speed and time, independently of
the production sampler.

Native execution requires a **clean regenerated CPU field-line build**, with
segment particle attachment enabled. Its headers and `build/libAMPS.a` must
match the source being tested; the fixture cannot repair an old mixed build.
Configure/build an owned checkout using the ordinary application setup (for
example `test/sep_parker_spiral__field_line-fte`). Set
`AMPS_FIELD_LINE_NATIVE_BUILD` to that build root. A configuration overlay must
explicitly declare the capability:

```json
{"schema_version":1,"capabilities":{"native-field-line":{"available":true,"policy":"required"}}}
```

Then run either registered native test, replacing the example paths:

```sh
AMPS_FIELD_LINE_NATIVE_BUILD=/path/to/clean/build-root \
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --config /path/to/native-capability.json \
  --test CORE-PIC-FIELD-LINE-NATIVE-SERIAL-001 --output-dir /tmp/amps-native-serial

AMPS_FIELD_LINE_NATIVE_BUILD=/path/to/clean/build-root \
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --config /path/to/native-capability.json \
  --test CORE-PIC-FIELD-LINE-NATIVE-MPI-001 --output-dir /tmp/amps-native-mpi
```

The serial test substitutes only the one-process MPI dependency operations and
prints that scope explicitly; unexpected communication aborts. It is **not MPI
qualification**. The other test starts real MPI with one rank, also not a
decomposition test. MPI startup failure is ERROR, not a physics PASS. Neither
replaces coupled SEP reference comparisons, multi-rank migration/decomposition,
or moving-geometry checkpoint validation. Default unavailable native capability
does not weaken existing gates; when declared required, a missing/mismatched
build or execution failure must fail the run.

The shared runner preserves transcripts, source/build hashes, native data and
the `.sampling` lifetime sidecar. Never regenerate a regression reference from
this manufactured fixture or count a dependency-adapter pass as coupled physics
validation.

## Existing interpolation, exchange and splitting scope

This suite qualifies the linear interpolation now used by production segment
getters, piecewise-curved arc interpolation, unique indexed payload exchange,
exact-intrinsic particle splitting, and conservative overlap remapping when a
line is repartitioned. Remote MPI transport remains a T10 gate; T8 checks the
packing/ownership and conservation algebra on one rank.

The split and invalid-input cases also exercise the exact clone kernel used by
`cFieldLineSegment::WeightedParticleSplitting`. They check relativistic kinetic
energy and unchanged phase space, opaque-history copying, the actual
parent/child callback arguments, equal-weight selection under permuted particle
handles, duplicate/zero identity rejection, subnormal-weight skipping, and an
application preflight rejection before allocation or weight changes. The
optional `FledLine::ParticleIdentity`, `PrepareExactSplit`, and `PostExactSplit`
hooks are registered before stepping. The key supplies the stable tie break;
preflight validates application lineage/state without consuming randomness;
post-split initializes declared application metadata after the exact clone.

These dependency-light records do not qualify PIC buffer storage, native MPI
transfer, or an application's child RNG/hazard policy. `CloneParticle` remains
byte-exact for migration and restart. Native application tests must separately
verify registration and complete stochastic lifecycle behavior. The field-line
weighted merger has a different conservation contract and is not qualified by
these split checks.
