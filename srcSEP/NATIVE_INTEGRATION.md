# D03 native AMPS integration and scientific execution gate

## Purpose and claim boundary

D01 and D02 can verify the SWCME adapter, canonical configuration schema, and
fingerprints without linking PIC. They cannot prove that the configured AMPS
application enters a production particle mover, maintains one SWCME generation
across ranks, or reproduces a checkpointed calculation. D03 supplies that
missing native evidence. It launches the linked `amps` executable; it never
reimplements a mover or substitutes a dependency-light driver.

The native C++ registry also exposes extended ID `D03PRE`. It is a fast linked
preflight used by `test/run_tests.py --all`: the callback requires exactly the
`parker`, `fte-dmumu`, and `fte-mfp` mover contracts and verifies two monotonic
SWCME prepared-state generations plus a canonical fingerprint. This catches a
mislinked executable before a site campaign is attempted. Its PASS is not D03
completion; it does not launch MPI ranks, advance a particle, write/reload a
checkpoint, or compare decomposition artifacts.

The development gate may report `SKIP` when a configured AMPS tree or reviewed
site campaign is unavailable. The release gate treats the same condition as
`INCOMPLETE` and returns nonzero. A source-only pass is therefore not presented
as native or scientific validation.

## Build contract

`test/run_native_integration.py` receives the enclosing AMPS source root and
its generated `Makefile.conf`. In release mode `--rebuild` is mandatory and
`--no-build` is forbidden. The runner first invokes the enclosing `make clean`,
then calls the srcSEP `strict-production` target. That target delegates to the
normal AMPS application build and audits `build/main/mainlib.a`; it does not
compile PIC-dependent translation units in the source application directory.

The evidence record stores the exact build commands, return codes, elapsed
times, complete logs and log checksums, compiler/MPI launcher versions, the
linked executable checksum, and the `Makefile.conf` checksum. `MAKEFLAGS`, such
as `-j16`, is inherited by both build commands.

## Native campaign manifest

Copy `test/native_integration_manifest.example.json` to a site-owned file,
remove `"template_only": true`, and replace every placeholder argument and
artifact path. The example is a schema and coverage template, not a runnable
AMPS input: `--site-input` and
`--checkpoint` deliberately mark places where the site must supply its actual
post-compile input/restart interface. Do not add `--particle-mover`, `--mover`,
or `--sep-mover` to case arguments; the gate reserves all aliases and appends
the canonical mover from the case record so the manifest cannot accidentally
test a different algorithm.

Each case declares:

- a filesystem-safe ID, kind, canonical mover, and positive MPI rank count;
- site arguments and optional environment values;
- expected process exit code;
- `minimum_source_state_id >= 2`, proving that the initial SWCME state was
  refreshed rather than merely constructed once;
- `minimum_particle_dispatches >= 1`, proving that a validated particle record
  reached and returned from the selected production mover; and
- required artifacts plus optional exact-equivalence groups.

`{amps_source}`, `{output_dir}`, `{case_dir}`, and `{executable}` placeholders
are expanded in argument and environment values for each case. `{ranks}` is
expanded in the MPI launcher. Environment variable names are literal.

The manifest validator requires, for each of `parker`, `fte-dmumu`, and
`fte-mfp`, one serial case and at least two distinct multi-rank decompositions.
One artifact group must span all three decompositions. It also requires
uninterrupted, checkpoint-writing, and resumed cases; uninterrupted and resumed
outputs must share an exact-equivalence group.

Choose equivalence artifacts that represent canonical physics state after
deterministic rank reduction. Do not compare timestamps, rank-local diagnostic
files, unordered particle dumps, or paths embedded in headers. Exact SHA-256
equality is intentional: D03 is checking decomposition/restart reproducibility,
not statistical similarity. If a scientifically stochastic product cannot be
serialized canonically, add a canonical reduced evidence product in the
application rather than weakening this gate to an undocumented tolerance.

## Runtime evidence

Every completed run must print:

- `RunConfiguration fingerprint=...`;
- `SWCME configuration fingerprint=...`;
- `final_source_state_id=...`;
- `completed_particle_dispatches=...`; and
- `mpi_consensus=pass`.

The dispatch counter is process-local during stepping and uses a relaxed atomic
because it is evidence, not physics state. It increments only after the
selected concrete mover returns. At shutdown the driver reduces it across MPI
ranks together with D01 query/recovery counters. A valid boundary deletion is a
completed dispatch; a validation or physics failure that exits before return is
not. Selection resets the counter, preventing an earlier mover choice from
supplying evidence for a later one.

The state ID and epoch are compared by MPI minimum/maximum before rank zero
prints the consensus marker. Run-configuration fingerprints must agree across
the decompositions of a mover, and one canonical SWCME fingerprint must cover
the entire campaign. Logs and every declared artifact are checksummed into
`native_integration_evidence.json`.

Use a new `--output-dir` for every executed campaign. The gate refuses a
nonempty `cases/` directory rather than deleting it or accepting a stale
artifact as current evidence.

## Commands and outcomes

From `AMPS/srcSEP`, a development probe is:

```sh
env MAKEFLAGS="-j16" test/run_tests.py \
  --suite d03-native-integration \
  --amps-source .. --make-config ../Makefile.conf \
  --output-dir test_output/d03
```

It returns a structured `SKIP` when the site manifest is absent. A release run
must use a reviewed, site-specific manifest and a clean build:

```sh
env MAKEFLAGS="-j16" test/run_tests.py \
  --suite d03-native-integration \
  --amps-source .. --make-config ../Makefile.conf \
  --native-manifest test/native_integration_manifest.site.json \
  --rebuild --release --output-dir test_output/d03-release
```

`PASS` means the clean configured build, all 12 or more native cases, MPI state
consensus, mover-dispatch minima, SWCME refresh, fingerprint checks, and exact
artifact groups passed. `FAIL` means an executed build/case or comparison
failed. `SKIP` is allowed only in development. `INCOMPLETE` is a release
prerequisite failure. `ERROR` denotes malformed invocation or infrastructure.

## Current limitations

D03 does not create a scientifically reviewed mesh, particle population,
restart file, or cluster launcher. Those are site inputs and must be reviewed
for adequate particles, timestep, output cadence, and SWCME validity interval.
The gate establishes execution and reproducibility evidence; it does not by
itself establish agreement with an observational event or an independent
transport solution.
