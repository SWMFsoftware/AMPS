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
floor.  The maintained RH solve brackets the nontrivial compressive root and
converges on bracket width.  It does not stop merely because the normalized
energy residual is small: close to `Mf=1` that residual is flat between the
identity and physical roots, so an early residual-only stop can satisfy the
flux tolerance while failing entropy or downstream-characteristic ordering.
The final reconstructed state must still pass the configured conservative
residual, compression, entropy and branch checks.  A root below the declared
minimum resolvable compression remains a typed unresolved weak shock.

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

### Triangular surface and physical patch measure

Surface quadrature is uniform in the generating-sphere normal coordinate
`mu=d.n`, whose area element is `a^2 dmu dphi`.  This avoids the singular
heliocentric-ray metric at SSE tangency.  The maintained
`triangular-sse-cap-v1` topology uses `polar_cells` intervals between
`mu_min=-sin(lambda)` and `mu=1`, and `azimuth_cells` periodic sectors.  It
stores every interval boundary through the last finite ring and stores the
`mu=1` apex exactly once.  Thus the front is a topological disk with one
physical boundary at the finite-width flank: its azimuth seam and apex are
closed, but no fictitious rear or side shock is added.

Every non-apex chart rectangle is split on a deterministic alternating
diagonal.  Its two triangles each own the exact curved area

```text
A_triangle = 0.5 a^2 Delta(mu) Delta(phi).
```

An apex fan triangle owns the full final chart-cell area because the complete
`mu=1` edge maps to the one apex point.  The exact areas sum, with compensated
accumulation, to `2 pi a^2 (1+sin(lambda))`.  The smaller planar chord area is
stored separately only to measure visualization-geometry convergence; it is
never used in an area ledger or physical normalization.

Coordinates are vertices, but shock state is face-centred.  Each triangle has
one stable ID, exact area, analytical representative position/normal/speed,
ambient query, status and optional RH state.  A super-fast/sub-fast boundary
is therefore never manufactured by averaging vertex flags.  This structure is
suitable for later deterministic injection weighting and surface sampling,
while the shared provider remains independent of particle allocation. The
srcSEP3D application now uses each patch's exact curved area for its physical
rate and square-root barycentric sampling on the corresponding maintained
planar triangle for native AMR ownership. Thus the injected point lies on the
discrete triangular shock representation, slightly inside the analytical SSE
surface. This is a declared geometric approximation whose magnitude is
measured by curved-versus-chord area convergence; it must not be described as
exact uniform sampling of the analytical curved patch.

Independent tests reconstruct edge incidence, orientation, connectivity,
Euler characteristic, analytical area closure and chord-area convergence.
Angular refinement and accepted-patch coverage still require convergence
evidence because a near-critical face can change physical classification.
The production positive-event regression also crosses a complete face ring at
`Mf=1.00106`: it requires a resolved weak jump before the crossing and a
physical `subfast-front` result afterward, with zero numerical-unknown area.

### Accepted incident-flux normalization

The shared diagnostics expose a side-effect-free normalization query at an
explicit requested apex radius. It first brackets the unique time on the
monotonically outward analytical trajectory, then calls the same production
surface/RH path used for committed epochs without replacing `Current()`. For
species population `s`, it evaluates

```text
Ndot_s = sum_{status=solved-fast-shock}
         n_{s,1} (V_n - U_1 dot n) A_curved .
```

The sum uses each face's exact curved area and actual upstream EOS density.
Sub-fast, non-forward and below-inner-boundary areas remain excluded physical
area rather than being renormalized. Any numerical-failure area rejects the
rate. Electron and proton rates are independent upstream populations; alpha
uses the event's declared alpha/proton number ratio. This is the gross number
of ambient particles swept through accepted front faces per second. It is not
an acceleration efficiency, seed-selection law or downstream volume source.

`srcSEP3D` combines this rate with its post-allocation global step and requested
model-particle count as `W_s=Ndot_s dt/N_model`. The shared provider itself
remains independent of particles and MPI. Test `RSH41` recomputes the face sum
independently from a production epoch, checks radius-domain rejection, and
proves that the query does not commit a generation.

For the currently implemented constant-weight application source, every rank
constructs the same keyed Poisson candidates at rate `Ndot_s/W_s`, but only the
rank owning the sampled active leaf allocates a particle. Face choice is
proportional to the per-face rate above. The momentum-number law is
`dN/dp proportional to p^(2-q)`, with either a fixed configured `q` or the
local test-particle DSA value `q=3X/(X-1)`, and the initial direction is
anti-sunward. These are srcSEP3D policies, not new provider physics; in
particular they do not add an efficiency model or a downstream plasma volume.
Because provider and immutable epoch cross the application boundary as
separate objects, the sampler requires the epoch event identity to equal the
provider's transitive physics fingerprint before it reads composition or
accumulates any face rate. A foreign or missing identity is a transactional
configuration conflict, not a zero-rate source.

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
- JSON epoch metadata includes vertex/triangle counts. CSV surface records
  include triangle connectivity, exact curved and planar chord areas, SI units,
  stable identifiers, absolute area ledgers, status and validity. CSV leaves
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

## Event file

A reduced shock-front event (`*.event`) is written in the generic AMPS runtime
input format (`src/INPUT_FILE.md`) as one
`#begin physics: shock_front` ... `#end physics` block.  `!` starts a comment
(there is no `#` comment) and `\` continues a line.  The keys are unchanged
from the former `key=value` record: the dotted vocabulary
(`run.start_s`, `ambient.reference_radius_m`, ...).  Inside the block a nested
block is the same key space as the dotted prefix, so

```text
#begin physics: shock_front
  #begin ambient
    reference_radius_m = 1.495978707e11   ! same as ambient.reference_radius_m
  #end ambient
#end physics
```

is equivalent to writing `ambient.reference_radius_m = 1.495978707e11`
directly in the block (`AMPS::Input::Flatten`); supplying the same key both
ways is a duplicate and fails.  The value shown is only a syntax illustration.

`ResolveConfiguration(bytes, reader, displayName)` parses that text with the
header-only reader `src/general/input/amps_input.h` (C++17, standard library
only, no PIC/MPI; `test/check_architecture.py` allows
`#include "amps_input.h"` as the one AMPS include).  Any other top-level block,
a `physics:` block other than `shock_front`, unknown keys, missing keys,
duplicate keys, unsupported enum values, non-SI values and a changed harmonic
checksum fail before the provider is constructed; diagnostics are prefixed
with `displayName:line:` (the file path, or a label for application-assembled
text).  `assets.harmonics_file` is a relative path whose bytes the caller
supplies through `reader` (the caller owns the filesystem policy); its
SHA-256 bytes enter the event fingerprint.  srcSEP and srcSEP3D resolve it
from the directory of the file that contains `#begin physics: shock_front`
(the event file itself when it is `#include`d into a srcSEP3D deck).  The physics identity is the
normalized manifest of the resolved assignments, so the change of file format,
comments and layout does not change it.

A stand-alone event file is used by the shared tests, viewers and srcSEP
`--reduced-shock-event`.  In `srcSEP3D` the event is `#include`d into the
application deck (`#include "handoff_smoke.event"`) instead of being named by
the retired `[background] model_asset` key, and the application passes the
canonical text of that block to the resolver.  An old `key=value` event can be
converted with `python3 tools/convert_amps_input.py event old.event -o
new.event`.  Runnable examples and the accepted application spelling are under
`srcSEP3D/examples/shock-front/`.

Schema `shock-front-ambient-v1.2` adds `history.model = fitted-quadratic`
to the existing `fixed-sse-pfss-parker-drag-v1` profile. It uses the fitted
HCI axis, fixed half-width, initial apex radius and speed, and
`history.apex_acceleration_m_s2`. The initial fit and dated ambient share
the UTC reference and start at zero elapsed seconds. `history.fit_end_s`
is the last admissible fit epoch. The matching drag branch starts with
the radius and speed calculated there, above
`history.minimum_handoff_radius_m`, and remains valid for
`history.drag_horizon_s`. A run that exceeds the fit plus drag coverage
is rejected. This schema derives both handoff radius and speed and omits
`history.final_pulse_speed_m_s`, `history.acceleration_duration_s` and
`handoff.apex_radius_m`. The SHA-256 provenance of the fit receipt is
recorded in `history.prepared_event_sha256` and included in the event
fingerprint. The native regression evaluates this provider on the same
physical clock used for particle transport.

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
- `numerics.*` owns surface resolution, the required
  `surface_topology=triangular-sse-cap-v1`, and weak/RH residual tolerances.
  The topology spelling enters the normalized event fingerprint because its
  facets and stable IDs are future source identities, not an output preference.
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
path.  The current registry contains 97 rows.  The retained unified native
receipt under
`test_output/reduced-front/runner/20261005T-native-restart-qualification-02/`
contains 93/0/0/0 because it predates the four rows added with `RSHAPP04`.
That boundary has separate passing evidence, but a new complete native run is
required before a 97-PASS aggregate can be claimed.
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
./Config.pl -application=sep3d   # deck selects runtime mover dispatch
make -j
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


## Read-only generating geometry for field-line consumers

`Provider::Sphere(epoch)` exposes the existing finite-SSE generating centre,
radius, unit direction and minimum supported ray cosine in HCI/SI. It calls
the maintained trajectory and does not prepare or publish an epoch. A sphere
root is only a candidate: consumers must use `EvaluateRay` to select the
leading supported surface and `EvaluateFrontPoint` for classification/normal/
speed/RH data. Rear-sphere roots never acquire physical shock support through
this API. This is an integration interface to the existing model, not a new
shock equation or downstream volume capability. The srcSEP SINT10–13 and
SINT13N gates exercise its finite support, tangency and material-time mapping.

## Schema v1.3: heliocentric sphere over a PFSS/Parker ambient with an open-tube wind

`schema = shock-front-ambient-v1.3`, `profile = heliocentric-sphere-pfss-parker-pchip-v1`.
It is used by the 2013-04-11 real-event test (`srcSEP/test/event_2013_04_11`),
whose README gives the derivations and the references.

**Front.** A closed Sun-centred sphere, `geometry.model = heliocentric-sphere`.
Its radius `R(t)` is the monotone PCHIP through `history.knot_times_s` and
`history.knot_radii_m` (`history.model = data-driven-pchip`, the SWCME
`data_driven` slopes). The normal is radial and `V_n = dR/dt`. The mesh is the
closed triangular sphere (`2 N_az (N_pol-1)` faces). No extrapolation is done:
the knots must start at `run.start_s` and reach `run.end_s`.

**Ambient.** WSO-convention PFSS harmonics below `ambient.source_surface_radius_m`
and the Parker continuation above it. The model has three additions:

1. **Wind.** `ambient.wind_speed_at_reference_m_s` selects the isothermal
   Parker (1958) transonic wind. Its sound speed `c_w` is the unique root of
   `U(r_ref; c_w) = U_ref`. The wind sets only `U` and `rho`; species
   temperatures, pressure and every Rankine–Hugoniot speed keep the declared
   `plasma.*` values.
2. **Magnetic calibration.** `ambient.magnetic_calibration_*` applies one
   magnetogram factor, so that `|B|` at a declared point and epoch equals a
   declared value. B is linear in the coefficients; topology, `U` and `rho`
   do not change.
3. **Open tubes.** `ambient.open_wind_density_model = isothermal-flux-tube` is
   the only value v1.3 accepts. In every open PFSS tube the plasma flows along
   B in the frame that corotates with the Sun:
   `U = Omega x r + s u b`, where `s` is the sector of the tube's foot.
   - *Induction.* The potential field is static in that frame, so the rigidly
     rotating field solves ideal induction exactly:
     `curl((Omega x r) x B) = -(Omega x r).grad B + Omega x B`.
   - *Continuity.* `rho u / |B| = K = Fm / (Rss^2 |B_r,ss|)`, the mass flux per
     unit magnetic flux at the tube's source-surface foot.
   - *Field-aligned momentum.* `u` is the subsonic root of
     `y - ln y = 4 ln(Rss/r_c) - 3 + 2 ln(|B_r,ss|/|B|) + 4 r_c/r`, with
     `y = u^2/c_w^2`. This is the isothermal Bernoulli integral along a tube of
     area `Phi/|B|` (Kopp & Holzer 1976). It reduces to Parker's relation where
     `|B| = |B_r,ss|(Rss/r)^2`.
   - `B`, `U` and `rho` are continuous with the exterior at `Rss`.
   - The centrifugal term is left out, as in the exterior wind. It is at most
     `Omega^2 Rss^3 / (2 GM)`, about 2e-4 of the gravitational potential.
   - Rejected, never clamped: a point with no subsonic root, `r_c <= Rss`, and
     the combination with `spherical-radial-mass-flux`.

**Evaluation.**
- The v1.3 ambient evaluates the magnetogram with
  `PfssHarmonics::EvaluateFieldCartesian`, which uses recurrences and agrees
  with `EvaluateCartesian` to about 3e-15.
- It finds an open tube's foot by integrating `dx/dr = B/B_r` in uniform radius
  steps of at most `ambient.trace_step_m`. This is valid while `B_r` keeps its
  sign and `|B_r| >= 0.2|B|`. Any other point uses the general arc-length
  `Trace`.
- `AmbientModel::SpeedBoundMPerS()` is a conservative bound on `|U|`. srcSEP
  uses it only to exclude line segments that cannot reach the front.
- Earlier profiles keep their historical evaluator bits.

**Material qualification.** srcSEP maps v1.3 events to the material profile
`reduced-pfss-open-parker-v1`
(`SEP::Background::Material::Profile::ReducedPfssOpenParkerV1`). It admits the
Parker exterior and `CoronalOpen`. Closed tubes, the sheath, ejecta and any
downstream volume remain unqualified.

**Tests.** `make shock-front-test` runs `test_shock_front_sphere`, HSP00–HSP12:
- HSP01–HSP05: PCHIP against SWCME `data_driven`, analytic derivatives,
  admission, exact sphere mesh and front membership.
- HSP06–HSP08: wind fit and decoupling, linear calibration, legacy closures.
- HSP09: along one traced open tube, `u' || B`, `rho u/|B|` and the Bernoulli
  constant are constant to about 2e-14, the flow is subsonic, and B/U/rho are
  continuous at `Rss` to about 1e-8.
- HSP10: the rejection cases.
- HSP11: recurrence against per-term PFSS evaluation.
- HSP12: the radius-parametrized foot against the general trace (`K` agrees to
  6.5e-9 with 2e7 m steps).
