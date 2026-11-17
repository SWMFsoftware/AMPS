# Provider-defined coefficient radius

Radius-dependent transport coefficients consume the same background authority
as B/rho/U and their gradients. They must not independently substitute the norm
of a global coordinate for a provider's declared radial coordinate.

`LocalBackgroundView` now supplies `heliocentricRadiusM`, an explicit known flag,
`radiusOriginM`, `heliocentricPositionM` and `radiusFrameIdentity`.
`EvaluateLocalBackgroundAt` establishes the frame before choosing a derivative
implementation:

| Provider | Radius |
|---|---|
| Analytic | `|x-origin|`, with origin from immutable initialization / the documented legacy declaration |
| SWCME | Existing global heliocentric frame, `|x|` |
| Reduced shock | Existing global heliocentric frame, `|x|` |
| SWMF/local imported handoff | Unavailable until that provider declares its radial frame |

A derivative fixture replaces grad(U), not radial origin. It neither clears a
known analytic radius nor adds endpoint-U validation. Conversely, explicit
grad(U) and finite imported coordinates do not establish a Sun/origin for a
radius-dependent coefficient. A radius that is unknown or undefined is not a
physical zero; `BuildSourceBoundInput` returns `UnsupportedConfiguration` before
publishing a coefficient. No global-origin fallback or inferred frame rotation
is introduced.

The optional parallel-diffusion library requires an already solar-relative
Cartesian position. Its native bridge now consumes `heliocentricPositionM`
from this same provider view instead of passing global native X unchanged.
The translation retains the Cartesian vector basis, so B and U are not rotated.
The Parker-facing spatial provider preflights a missing frame as typed
`UnsupportedConfiguration` before entering its numerical derivative wrapper;
that capability error must not be hidden as a generic failed stencil. No
canonical library formula or optional isotope/decomposition input is changed.

The single immutable semantic identity is
`Background::RadiusFrame::ContractIdentity` in
`util/sep_provider_radius_frame.h`. The existing background fingerprint owner
includes it; analytic origin/speed values still come from the original
initialization owner. Coefficient provenance carries the actual radius [m] and
frame identity. The native particle / physics-checkpoint wire layouts do not
change, but an older physical interpretation has a different configuration
fingerprint and cannot silently resume as the corrected model.

`RADIUSN01` is an Extended native callback. Its explicit
`native-startup-input-schema6` descriptor capability opts its individually
selected `--test-input` into the existing `Initialization::LoadFile/Install`
startup pipeline. Production `--input` still rejects test mode; other tests'
input manifests are never guessed to be startup records. There is no mutable
initialization override or partially filled dummy configuration. The dedicated
runner derives a complete owned schema6 record
from `examples/sep_parker_mesh.in`, explicitly declaring an offset origin,
725 km/s wind, 17 points, a 3e11 m outer radius and a coarse unrefined support
mesh. It retains the reviewed source, observer and thermal fields; existing
`--shock-model analytical` is the explicit diagnostic provider override, and
`--particle-mover parker` retains explicit CLI precedence through file layering. The
unchanged 2e11 m arclength is strictly bounded by the outer domain. This is an
initialization/read-only accessor test, not transport or mesh convergence.

The callback reads the actual native position/B/U, independently derives
`r=|x-origin|` and `divU=2U/r`, and calls the actual QLT1 mean-free-path provider.
Its independent reference is `lambda=Lc/(deltaB/B)^2*(rL/Lc)^(1/3)`, with
`Lc` proportional to the declared radius, actual compiled species mass/charge,
and explicit controlled coefficient parameters. The old global-zero radius
must produce a materially different negative control. Post-fix checks also
cover direct background radius/frame, provenance, the derivative-fixture branch
and explicit rejection of an unknown SWMF origin despite a supplied tensor.
Coefficient configuration, immutable startup, snapshot identity and particle
ownership are preserved; no source or host step is run.

`RADIUSN02`, in the same factory/group, separately installs the actual optional
library through complete immutable startup. Its reviewed base is
`examples/sep_parker_mesh_parallel_diffusion.in`; its explicit diagnostic law is
`lambda=3e9 m * (v/(2e6 m/s))^0.75 * (r/(1 AU))^1.5` with
`kappa=v*lambda/3`. The callback independently evaluates this law at the native
solar-relative position and requires a materially different global-zero
negative control. It checks actual library selection/configuration, finite
refined derivative, unknown imported-origin rejection and unchanged native,
initialization, library and background ownership. This is a software/coordinate
contract fixture, not a fitted or observationally validated transport law.

Prepare a reviewable input without claiming a test pass:

```sh
python3 -B srcSEP/test/run_native_coefficient_radius_tests.py \
  --amps ./amps --output /tmp/sep-radius-review-new --prepare-only
```

Execute the separate process fixture:

```sh
python3 -B srcSEP/test/run_native_coefficient_radius_tests.py \
  --amps ./amps --mpi-ranks 2 --output /tmp/sep-radius-native-new
python3 -B srcSEP/test/run_native_coefficient_radius_tests.py \
  --amps ./amps --case RADIUSN02 --mpi-ranks 2 --output /tmp/sep-library-radius-native-new
```

Output directories must be new/empty. Input/template/executable hashes, commands,
native logs, structured reports and initialization products are retained in
owned scratch. The input's relative output leaves are written only in its
owned native working directory; no initialization-only output flag is used
in test mode. SKIP is not a pass. The isolated BUILD10a old-producer build
reproduced both intended failures in actual callbacks on two MPI ranks, with
job exit 1 and no startup/harness errors. The expected provider radius was
`2.0162684853514194e10 m`, while the global-zero norm was
`4.153437523037244e10 m`. QLT1 returned the old-global result
`5.043862336867351e9 m` rather than the independent
`3.1154707493976245e9 m`. The optional library returned
`8.698631001278114e13 m²/s` rather than `2.9421304733960125e13 m²/s` and failed
the typed unknown-origin control. Earlier BUILD10 startup rejection did not
reach the callback and is retained separately; it is not the physics failure
evidence. Post-fix native qualification remains pending.

Source syntax checks and runner preparation do not establish the post-fix
gate. The runner's existing default case
and summary schema/keys remain stable; `--case RADIUSN02` selects the second
explicit input protocol. Source CLI tests pass 10/10 and parallel-diffusion
binding tests pass 5/5; actual LoadFile/Configure accepts the owned library
input. Those checks do not establish native coordinate qualification. Track
actual evidence in STATUS.

No isotope/nucleon metadata guess, imported composition/temperature law,
material-grid update, kinetic source matching or full production closure is
resolved by this radius repair.
