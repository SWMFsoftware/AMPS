# Background providers and immutable samples

This directory converts model state into the fields consumed by srcSEP3D.
Providers and snapshot construction have no PIC or MPI dependency. The AMPS
adapter in [`../main_lib.cpp`](../main_lib.cpp) owns center-node storage,
collective validation, ghost exchange and publication. Model adapters must not
write a cell or call an MPI collective from `Evaluate`.

## Source responsibilities

| File | Responsibility |
|---|---|
| `bg_provider.h/.cpp` | Interface, metadata, capabilities, per-point batch fallback and vector-derivative algebra |
| `bg_parker.h/.cpp` | Canonical SWCME Parker/Leblanc ambient state with analytic derivatives |
| `bg_swcme.h/.cpp` | Canonical 3-D regional SWCME primitives and numerical vector derivatives |
| `bg_swmf.h/.cpp` | Validated host-imported state; acquisition remains owned by the coupled host |
| `background_snapshot.h/.cpp` | Complete-sample validation, immutable candidate construction and optional bracketed snapshot interpolation |

[`../BACKGROUND_FIELD.md`](../BACKGROUND_FIELD.md) describes the full storage
and field equations. [`runtime/README.md`](../runtime/README.md) describes how
a provider is selected and published. The
[SWCME update guide](../../SWCME_MESH_BACKGROUND_CHANGES.md) gives build and run
commands and the boundaries of the implemented physics.

## Preparation, evaluation and ownership

Construct a provider from frozen parameters, call `Validate`, then call
`Prepare(timeS)` once at an update boundary. Preparation must construct a
complete candidate before replacing its previous state. On success,
`PreparedMetadata` exposes the epoch, validity interval, frame, provenance,
configuration fingerprint and positive generation. Its non-owning pointer may
be invalidated by the next successful preparation; snapshots copy metadata.

Reuse `Evaluate` or `EvaluateBatchDetailed` at that prepared epoch. Evaluation
must not advance time, run a process, acquire new source data, or update
metadata. Do not overlap `Prepare` with evaluation on the same instance.
Successful batch points replace their output; rejected points receive an
individual status and preserve their previous output. The aggregate status is
the first failure. A provider-specific override must retain these semantics.

`BackgroundSnapshotBuilder` evaluates into temporary storage and assigns the
output shared pointer only after every point succeeds. Previously published
snapshots stay immutable when a provider advances. An empty local position
array is legal: a pruned MPI rank still carries the prepared identity and joins
publication. The AMPS boundary requires positive **global** physical coverage.

## Units and derivative convention

| Sample value | Meaning and units |
|---|---|
| `B`, `absB`, `bHat` | Vector IMF [T], magnitude [T], unit direction |
| `U` | Plasma bulk velocity [m/s] |
| `gradB(i,j)`, `gradU(i,j)` | Partial derivative of component `i` with respect to coordinate `j` [T/m, 1/s] |
| `numberDensityM3` | Electron density [m^-3] in the current canonical closures |
| `temperatureK`, `pressurePa` | Proton temperature [K], total thermal pressure [Pa] |
| `alfvenSpeedMpS` | `absB/sqrt(mu0*rho)` [m/s], with configured composition |
| `divU`, `fieldAlignedStrain` | `trace(gradU)` and `bHat*bHat:gradU` [1/s] |
| `divBhat`, `curvature` | Divergence of the unit field and `(bHat dot grad)bHat` [1/m] |
| `focusingLenM` | `-1/(bHat dot grad(log(absB)))` [m] |

`CompleteVectorDerivatives` needs a non-null sample, finite nonzero `B`, and
complete `gradB`/`gradU`. It derives the transport scalars and vectors but does
not validate the result or set its status, generation or valid flag. In
particular, it uses

\[
\nabla\cdot\hat b=\frac{\nabla\cdot B}{|B|}
 -\hat b\cdot\nabla\log|B|,
\]

without assuming that a prescribed regional field is solenoidal. Zero
directional magnitude derivative gives positive infinite focusing length;
the completeness validator accepts that limit only with a finite tensor and
an exactly zero sampled derivative. Other infinite primitives are rejected.

New availability flags (`hasGradB`, `hasGradU`, `hasDivBhat`, `hasCurvature`)
declare supplied derivatives independently of analytic provenance. Numerical
or imported derivatives must pass the same finite checks. Existing analytic
flags continue to imply availability.

## SWCME adapter choices

`SwcmeBackgroundProvider` uses the canonical checked primitive API, including
RH-heated total pressure. It preserves configured electron/proton and
alpha/proton temperature ratios while scaling all temperatures to recover
that pressure. This is a declared heating partition, not an independently
solved species energy equation. SHOCK_ONLY samples ambient fields; FULL_ICME
samples the finite shock, sheath and simple ejecta prescribed by SWCME.

The 1..1.05-Rs physical shell in the example uses a matching ambient
continuation below the canonical CME domain. `Prepare` rejects overlap of the
entire CME trailing transition with the handoff. This continuation does not
provide a coronal solar-wind model or a corona-to-SWCME coupling interface.

Cartesian second-order derivatives use seven points per cell. A relative
step of `1e-4*r`, initially floored at 1 m, is capped by shock width/16.
At the excluded solar surface, an outward second-order stencil replaces a
central stencil that would enter the Sun. Scratch batches contain at most
256 cells. On an exceptional canonical batch failure, individual stencils
are rechecked to retain valid points and precise statuses. The complete mesh
candidate still fails if any required cell fails.

Canonical failures retain their full status reason, batch index, epoch and
heliocentric position in metres. This includes stencil positions, so a rejected
mesh initialization can be reproduced without guessing the owner-cell centre.
The canonical fast-branch/switch-on solver revision also enters the manifest.

The relative step, layer cap, heating partition, inner continuation and
coordinate origin enter the resolved manifest. They are part of physics
identity, not undocumented numerical defaults.

## Tests

From the AMPS root:

```sh
python3 srcSEP3D/test/run_tests.py --suite phase-b --rebuild --output-dir test_output/background
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in --test-steps 2
```

`SWBG3D01–08` exercise real canonical providers without AMPS/MPI, including
polar region snapshots with all derivative stencils and an independent analytic
parallel shock limit. See [the polar initialization fix](../../SWCME_POLAR_SHOCK_FIX.md).
The second
command discovers `SWBGAMPS01–03` from the linked executable and checks owner
buffers, actual refresh counts and representatives in received remote blocks.
Their prerequisite SKIPs are documented in
[`../test/README.md`](../test/README.md). Passing publication checks does not
establish observational skill or resolved particle-acceleration convergence.

## Finite-SSE support

`SwcmeBackgroundProvider` accepts Sphere and SSE. The canonical model returns
ambient Parker/Leblanc primitives outside the finite cap; the adapter does not
clip directions into an artificial flank. Regions and derivative step caps use
the **local** shock radius, so a narrow flank layer does not inherit an apex-sized
finite-difference step. Preparation checks the minimum radius at the tangent
flank, including the trailing transition, against the 1.05-Rs ambient handoff.
A failed epoch leaves the prior prepared state intact.

The finite cap boundary can terminate phenomenological FULL_ICME ejecta abruptly.
Finite derivatives near that edge depend on the stencil; this implementation
does not claim a smooth lateral MHD solution or divergence-free flux rope.
Use converged regions away from that truncation for acceleration studies.

Run `python3 srcSEP3D/test/run_tests.py --group SSE3D --rebuild` for the parser,
canonical geometry, background, source, crossing, subcycling and checkpoint prerequisites.
For native publication use the coupled runner with
`--test-input srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in`.


### Finite-SSE support and weak surface jumps (2026-10-02)

Canonical 3-D field evaluation now classifies local geometric support before
solving the surface RH state for each ray. Upstream, post-ICME and outside-cap
points are ambient by the prescribed model; all smoothing layers stay on the
regional path. This applies consistently to the n/V, vector and full primitive
interfaces. It prevents a far-upstream finite-SSE stencil from failing because
a surface shock on that ray lies near Mach one.

The stable weak-jump solver resolves the reported `t=180 s` state without a
Mach cutoff, angular masking or ambient substitution at the front. A genuine
in-CME solver failure still rejects the background candidate collectively.
The provider manifest now includes `weak-delta-v2` and the geometry-first
regional-support policy, so snapshots produced with different numerical
policies cannot silently share an identity. See
[the failure analysis](../SSE_WEAK_SHOCK_FIX.md).
