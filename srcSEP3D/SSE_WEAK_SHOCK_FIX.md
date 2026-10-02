# Finite-SSE runtime weak-shock correction — 2026-10-02

## Cause of the reported abort

The failure is reproduced with the supplied complete SSE input at `t=180 s`
and both reported positions:

```
(78429796951.338379, -31283096398.984375, -2999748969.765625) m
(78429796951.338379, -31283096398.984375,  2999748969.765625) m
```

The sampled radius is `8.4491796430e10 m`; the actual front along that ray is
only `1.2461375682e10 m`. The sampled mesh point is far upstream and must
have the prescribed ambient state. Its ray's surface fast Mach number evolves
from `0.9963560037` at 120 s to **`1.000000867696212`** at 180 s. The former
solver rejected every positive excess at or below `1e-6`, even when a stable
physical weak root existed. The 3-D vector field kernel's geometric ambient
shortcut applied only to Sphere, so an irrelevant surface solve was attempted
for this SSE ambient point and its derivative stencils. That combination
produced `NUMERICALLY_UNRESOLVED_WEAK_SHOCK` and rejected the rank-local
background candidate. The collective error then reached `StopWithStatus`
and intentionally aborted all ranks.

## Corrections

`src/models/swcme/swcme_shock.hpp` now routes representable weak shocks to the
analytically deflated fast-interval cubic, evaluated about `delta=r-1` with
long-double intermediates. The old divided energy-flux scan is bypassed for
`M_fast-1 <= 1e-4`; that value selects an algorithm and imposes no physical
Mach or compression floor. The explicit unresolved guard now covers only a
binary64 uncertainty margin of `128*epsilon(double)`, about `2.84e-14`.

The reported physical surface state converges to compression
**`1.0000011572267356`**. Its independently reconstructed mass, magnetic,
electric, momentum and energy fluxes pass the existing acceptance gates, as
do entropy and downstream characteristic checks. No empirical compression or
ambient replacement is used at the shock. A frozen 80-digit direct conserved-
flux solve gives `1.000001157226735592412583175971707701369547676320658815...`;
the audit-only verifier is `test/reference/verify_sse_weak_flank.py`.

`src/models/swcme/swcme3d.cpp/.hpp` now apply one geometry-only support test to
all three Cartesian field interfaces and all canonical shapes. The local
ray front, rather than SSE apex radius, defines shock and trailing transition
bounds. Ambient points skip the irrelevant RH solve. Every finite smoothing
layer remains on the full regional path. Direct shock/source diagnostics
still solve the actual shock, and genuine unresolved in-CME field queries
remain explicit failures with untouched rejected outputs.

`srcSEP3D/background/bg_swcme.cpp` versions both numerical policies in its
manifest. Detailed comments, canonical/application READMEs, and the complete
SSE example document the behavior. The input's physical parameters are
unchanged; new comments request the ten-step regression run.

## Rebuild and native rerun

The archive is a cumulative AMPS-root source overlay, including finite SSE,
the earlier polar jump/runtime DATAFILE corrections, and this weak-shock fix.
It contains directly modified files, with no patch installer or compiled
binaries. After extracting it, refresh any configured/generated application
source copies through the usual AMPS configure procedure. Force the canonical
shared SWCME archive to rebuild, then rebuild the application:

```sh
cd /nobackupp17/vtenishe/Mars1/AMPS
tar -xzf /path/to/SEP3D_finite_SSE_weak_shock_fix_20261002.tar.gz
env MAKEFLAGS="-j16" make -B -C src/models/swcme
env MAKEFLAGS="-j16" make amps
```

Run the native test into a fresh directory:

```sh
mkdir -p test_output/sse-weak-shock/native
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in --test-steps 10 --expect-mpi-ranks 4 --test-json test_output/sse-weak-shock/native/native.json --artifact-directory test_output/sse-weak-shock/native/artifacts
```

Portable application checks:

```sh
python3 srcSEP3D/test/run_tests.py --group SSE3D --rebuild --output-dir test_output/sse-weak-portable
python3 srcSEP3D/test/run_tests.py --suite standalone --no-build --output-dir test_output/sse-weak-standalone
python3 srcSEP3D/test/reference/verify_sse_weak_flank.py
```

## Verification scope

See [SSE_VERIFICATION.md](SSE_VERIFICATION.md) and evidence under
`test/evidence/weak_sse_20261002/`. SSE3D08 executes the exact mirrored
coordinates, neighboring points, full derivative stencils and the actual
shock layer at every 60-s epoch through ten updates. SSE3D09 verifies that
ambient support can succeed while deliberately roundoff-unresolved in-CME
queries and direct diagnostics still fail. SHK12/SHK15 extend numerical
coverage across obliquity, beta, gamma and 25,000 weak random states.

The target NASA four-rank HPE MPT executable has not been built or run here.
The successful portable checks reproduce and correct the numerical failure;
the native rerun above is still needed to verify the complete MPI lifecycle.
SSE's prescribed geometry and phenomenological interior retain the physical
limits documented in [SSE.md](SSE.md).
