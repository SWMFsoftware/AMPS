# Reduced shock-front plus ambient provider

This directory implements the version-1.1 reduced model specified by
`../docs/shock-front/model.md`.  It is independent of PIC, MPI, SEP sources and
particle state.  The production role is deliberately narrower than BG3D-4:
the provider returns an undisturbed ambient reference throughout its declared
domain and immediate one-sided RH limits only on accepted front samples.  It
does not produce a sheath, contact, ejecta, downstream volume or retained
shocked material.

## Implemented first profile

- SI values in an HCI Cartesian frame, outward normal from geometric interior
  toward upstream, and positive inflow `w1=Vn-U1.n`.
- Checksummed real-orthonormal PFSS harmonic asset, the maintained composite
  PFSS/Parker ambient implementation, and a fully specified multi-species
  ideal-gas EOS.  The first reduced profile explicitly selects spherical
  `rho U r^2` open-wind density so the section-18.5 polar reference is
  reproducible; the legacy field-aligned `rho U/B` option remains available
  and unchanged for the full-CME profile.
- Fixed-direction, fixed-width finite SSE leading surface.  For apex distance
  `Ra`, `c=Ra/(1+sin(lambda))`, `a=c sin(lambda)`, position is on the generating
  sphere, and local normal speed is
  `Vn=Va(d.n+sin(lambda))/(1+sin(lambda))`.  Consequently a fixed-width tangent
  edge has zero normal speed; it is geometric support, not accepted shock
  coverage.
- Constant-speed smoke history or the analytical quintic velocity pulse,
  followed by the maintained-compatible quadratic-drag analytical outer law.
  Radius and speed are matched at the bracketed handoff.  The fixture is C1;
  its one-sided acceleration can jump.
- Orientation-dependent fast speed and the maintained ideal-MHD RH solver.
  Non-forward, sub-fast, unresolved-weak and invalid-jump records retain no
  downstream primitive.  `criticality_model=none` is the only accepted first
  profile.
- Transactional epochs: failed required coverage never replaces the committed
  generation.  Sub-fast/non-forward classifications are valid science results
  under `report-only` and do not erase older records.

The production local-jump entry point evaluates
`w1=Vn-U1.n`, `Mf=w1/cf`, and only calls the maintained conservative solver
for `w1>0` and `Mf>1`.  Its downstream state satisfies the normal mass and
magnetic-flux jumps, tangential electric and momentum jumps, normal momentum,
and total-energy flux.  A positive Mach excess whose nontrivial compression
root cannot be separated from `X=1` is `numerically-unresolved-weak-shock`, not
`subfast-front`.  The first profile uses no Mach, speed, field, or compression
floor.

For an accepted jump the additional diagnostics are

```text
CB = |B2|/|B1|
vHT,t = u1,t - (u1,n/B1,n) B1,t
u1,HT = u1 - vHT,t = (u1,n/B1,n) B1
```

where `u=U-Vn*n` is shock-frame velocity.  `vHT,t` is a Galilean boost, not
the incident-flow speed `|u1,HT|=w1/|cos(thetaBn)|`.  Exact perpendicularity
has no finite conventional HT boost; a magnetic null has neither magnetic
compression nor direction.  Both remain valid hydrodynamic/MHD jump limits
with an explicitly absent diagnostic.  Near perpendicularity reports the
small normal-field cosine and large incident speed without an obliquity floor.

Surface quadrature is uniform in the generating-sphere normal coordinate
`mu=d.n`, whose area element is `a^2 dmu dphi`.  This avoids the singular
heliocentric-ray metric at SSE tangency.  Angular refinement and accepted-patch
coverage still require convergence evidence.

## Connectivity, observers, outputs, and restart

`diagnostics.h` operates on the same analytical provider and immutable epoch:

- Polyline/front intersections solve each segment's generating-sphere
  quadratic, collapse only a roundoff-size double-root discriminant, retain
  tangent roots, deduplicate shared vertices, and then re-query the leading
  finite SSE surface.  Thus a rear-sphere root cannot become a cobpoint.
  Attached polarity and local jump data come from the canonical ambient and
  production classifier.  The supplied polyline is an instantaneous ambient
  field-line approximation, not a particle orbit.
- Fixed or constant-velocity observer passages use bracketed signed-distance
  roots plus minimization of local `|g|` minima for grazes.  The published
  relative speed is `Vn-vobs.n`.  A geometric hit retains its local status, so
  an accepted shock hit, sub-fast hit, numerical unknown, graze, and finite-
  support miss are distinguishable.
- JSON epoch metadata and CSV surface records include SI units in column names,
  stable identifiers, absolute area ledgers, status and validity.  CSV leaves
  absent downstream values empty; it never serializes them as physical zeros.
- A background-only restart binds the event/asset fingerprint, committed clock
  and generation, phase, and exact handoff state.  Restore re-evaluates the
  analytical epoch and rejects a changed physical identity transactionally.
  That shared serialization is a physics-kernel reference only.  Native
  qualification additionally writes and reloads the AMPS application
  checkpoint, reconstructs the mesh-owned provider state, exchanges received
  ghosts and advances the normal AMPS time loop.
- Connected derivatives require the same stable branch on three ordered
  epochs.  They differentiate values at the moving root, including
  `xdot.grad(F)` implicitly; mergers, grazes, selection changes, or invalid
  samples return absence.  Ensemble accounting requires one common label and
  normalized declared weights, retains known no-shock and numerical-unknown
  weight, and reports probability bounds instead of renormalizing survivors.

## Event asset

The resolver accepts a strict, duplicate-free `key=value` file.  Unknown keys,
missing keys, unsupported enum values, non-SI values and a changed harmonic
checksum fail before the provider is constructed.  `assets.harmonics_file` is
resolved by the caller relative to the event file; its SHA-256 bytes enter the
event fingerprint.  Runnable examples and the accepted application spelling
are under `srcSEP3D/examples/shock-front/`.

The input groups are deliberately explicit:

- `run.*` owns the UTC reference, elapsed-SI support, HCI frame, background
  cadence and the hard `particle_mode=disabled` capability.
- `domain.*`, `ambient.*` and `plasma.*` own radial support, the checksummed
  PFSS/Parker state, density law, temperatures, composition and RH adiabatic
  index.  The ambient is queried at the cell's current position; it is not a
  quiet-state placeholder that changes behind the front.
- `geometry.*` owns the finite SSE direction and half width.  The surface has
  no artificial rear or side shock.
- `history.*` owns the prescribed low-coronal apex motion.  `handoff.*` owns
  the continuous analytical outer trajectory and the policy for scientifically
  valid sub-fast/non-forward states.
- `numerics.*` owns surface quadrature and weak/RH residual tolerances.
  `endpoint.*` owns the geometric target and observer, without implying that
  a shock must survive there.

See the commented `.event` and `.in` files beside the examples for concrete
values.  Changing any physical value or referenced asset changes the frozen
event/application identity.

## Application publication, epochs, and MPI ownership

The `srcSEP3D` adapter is intentionally thin.  It prepares one immutable shared
front epoch and exposes the maintained ambient primitive and Cartesian
derivatives through the existing `BackgroundProvider` ABI.  Cell classification
does not paint downstream values: every native cell receives the undisturbed
ambient reference, while front geometry and one-sided jump records remain a
separate diagnostic object.

For a prepared time `t`, the adapter derives a generation from the event's
declared background cadence, prepares the shared provider transactionally, and
only then publishes matching metadata and samples.  AMPS writes owner-cell
current/previous DATAFILE slots before its maintained halo exchange.  Received
physical ghost cells are independently read back against the provider.  A
rank-local candidate rejection is reduced collectively; any rejecting rank
leaves the committed epoch, generation, owner storage and front inventory
unchanged.  Event identity, front generation, ambient generation, apex state
and absolute area ledgers enter the MPI fingerprint.

The native zero-particle check reads the actual AMPS linked lists and source
ledger after stepping.  It requires `particles_per_cell=0`, source disabled,
global allocated particle count zero and global injected count zero; an empty
portable ledger alone is not evidence.

## Build and tests

The complete selected reduced-profile campaign has one srcSEP3D orchestrator:

```bash
cd /home/vtenishe/Mars2/AMPS
python3 srcSEP3D/test/run_reduced_shock_front.py
```

It runs the shared RSH/architecture gates, portable RSHAPP boundary, native
one-/four-rank smoke cases, actual four-to-four and one-to-four checkpoint/
resume comparisons, and the separate four-rank 1-AU campaign.  Fresh phase
logs, native receipts and `summary.{txt,json}` are written below
`test_output/reduced-front/runner/`; every FAIL/ERROR is printed with its log
path.  The current validated aggregate is 93/0/0/0, with evidence under
`test_output/reduced-front/runner/20261005T-native-restart-qualification-02/`.
On the 2026-10-05 host,
the 340-step native phase alone took about 29 minutes, so this command is not a
routine short test.  Use `--list` or `--dry-run` to inspect coverage/commands,
and `--skip-native` only when explicit native SKIPs are acceptable.

The runner derives `restart-smoke.in` from the handoff fixture, enables a
five-tick checkpoint cadence, and gives every process group a unique working
directory.  A resumed segment necessarily names the checkpoint as an input
and publishes into a different directory.  `RunConfiguration3D` therefore
stores two related identities: the full resolved manifest for provenance and
a restart-compatibility manifest containing all resolved physics/assets plus
the output/checkpoint cadence clocks, but no relocatable filesystem path.  The
physics fingerprint and native storage-layout fingerprint remain separate
mandatory comparisons.  Removing paths from the restart comparison does not
permit a changed event, model parameter, mesh layout or cadence.

At final tick 10, the native comparator requires exact reduced-front/event
identity, geometry, speed, per-record shock classification, epochs and
generations.  It independently fingerprints the actual owner-cell ambient
plasma/IMF values and checks received physical ghosts on four ranks, rollback,
checkpoint sequence, zero allocation and zero injected/remaining particles.
The four-to-four case qualifies same-rank restart; the one-to-four case
qualifies deterministic repartition for this background-only zero-particle
profile.  Nonempty particle repartition remains deliberately unclaimed.

Shared development build (not native AMPS evidence):

```bash
make -C src/models/sep_corona_swcme -j16 shock-front-test
src/models/sep_corona_swcme/build/test_shock_front
```

Native rebuilds must use the repository-root preflight recorded in
`CODEX_REDUCED_SHOCK_PLAN.md`; a source-only compile does not qualify owner or
received-ghost storage.  The maintained native workflow is:

```bash
cd /home/vtenishe/Mars2/AMPS
pwd
ps -C make -C gmake -C g++ -C gcc -C cc1plus -C mpiexec -C mpirun -C amps \
  -o pid=,ppid=,stat=,etime=,args=
rm -rf -- build
./Config.pl -application=sep3d
./ampsConfig.pl -input sep3d.input -no-compile
make -C srcSEP3D prepare-production
make -j16 amps
```

The process check must be empty before removal.  Do not remove model build
directories, site configuration, libraries or evidence.

Use a unique output root and override the deck's runtime product directory;
`--test-json` and `--artifact-directory` alone do not redirect observer output:

```bash
REDUCED_OUTPUT_ROOT=test_output/reduced-front/runs/UNIQUE_TAG
mkdir -p "$REDUCED_OUTPUT_ROOT/handoff-1/artifacts" \
  "$REDUCED_OUTPUT_ROOT/handoff-4/artifacts" \
  "$REDUCED_OUTPUT_ROOT/long-4/artifacts"

mpiexec -n 1 ./amps --test-suite sep-corona \
  --test-input srcSEP3D/examples/shock-front/handoff_smoke.in \
  --test-steps 10 --expect-mpi-ranks 1 \
  --output-dir "$REDUCED_OUTPUT_ROOT/handoff-1/products" \
  --test-json "$REDUCED_OUTPUT_ROOT/handoff-1/native.json" \
  --artifact-directory "$REDUCED_OUTPUT_ROOT/handoff-1/artifacts"

mpiexec -n 4 ./amps --test-suite sep-corona \
  --test-input srcSEP3D/examples/shock-front/handoff_smoke.in \
  --test-steps 10 --expect-mpi-ranks 4 \
  --output-dir "$REDUCED_OUTPUT_ROOT/handoff-4/products" \
  --test-json "$REDUCED_OUTPUT_ROOT/handoff-4/native.json" \
  --artifact-directory "$REDUCED_OUTPUT_ROOT/handoff-4/artifacts"

mpiexec -n 4 ./amps \
  --test RSH24 --test RSH25 --test RSH26 --test RSH27 --test RSH28 \
  --test-input srcSEP3D/examples/shock-front/corona_to_1au.in \
  --test-steps 340 --expect-mpi-ranks 4 \
  --output-dir "$REDUCED_OUTPUT_ROOT/long-4/products" \
  --test-json "$REDUCED_OUTPUT_ROOT/long-4/native.json" \
  --artifact-directory "$REDUCED_OUTPUT_ROOT/long-4/artifacts"
```

The smoke host step is 60 s, so ten actual steps include initialization, a
pre-handoff epoch, the exact 69.57 s crossing between ticks 1 and 2, and later
epochs.  The long host step is 600 s; 340 actual steps commit 204000 s, beyond
the independent continuous-time 1-AU root.

Address/undefined sanitizer verification uses an isolated temporary build:

```bash
make -C src/models/sep_corona_swcme -j16 \
  BUILD=/tmp/amps-reduced-shock-sanitized \
  CXXFLAGS='-O1 -g -std=c++17 -Wall -Wextra -Wpedantic -Werror -fsanitize=address,undefined -fno-omit-frame-pointer' \
  /tmp/amps-reduced-shock-sanitized/test_shock_front
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  /tmp/amps-reduced-shock-sanitized/test_shock_front
```

LeakSanitizer is a separate capability; it cannot run in a ptrace-managed
environment and must be reported unavailable rather than counted as PASS.

## Current validation boundary

The selected fixed-SSE/PFSS/Parker/drag profile is implemented and validated
through reduced stages S0--S7.  The 2026-10-04 evidence has shared 54/0/0/0,
srcSEP3D portable 146/0/0/0, one- and four-rank reduced smoke 4/0/0/0 each,
and the long native reduced selection 5/0/0/0.  The long trajectory reaches
1 AU geometrically at `203884.49378697205 s`; its exact polar observer result
is `non-forward-inflow`, so accepted shock arrival is false.  That is a valid
negative physical classification, not an accepted-shock success.

These results qualify only the declared reduced capabilities.  The provider
still has no sheath, contact, ejecta, shocked downstream volume, front dynamics,
particle source or event-specific observational validation.  Optional
driver/standoff, GSD and nonlinear-formation extensions RSH33--RSH35 remain
unimplemented.  Passing the reduced suite never qualifies the preserved,
deferred BG3D-4 full-volume model.  Exact hashes, failed attempts, suite SKIPs
and evidence paths are maintained in `CODEX_REDUCED_SHOCK_PLAN.md`.
