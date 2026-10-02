# Finite self-similar-expansion fronts in SEP3D

SEP3D accepts the canonical SWCME `sphere` and `sse` geometries. SSE supplies
an outward, finite angular cap; the spherical example remains a global
Sun-centered verification control. Ellipsoid particle crossing is still
unsupported and is rejected before AMPS allocation.

## Geometry and units

Set the following canonical assignments in a complete schema-3/4 input:

```ini
[swcme]
geometry.shape = sse
geometry.half_width_rad = 0.6981317007977318
geometry.cme_direction_x = 1
geometry.cme_direction_y = 0
geometry.cme_direction_z = 0
```

The example is a 40-degree half width, or 80-degree total opening, propagating
along the simulation frame's +X axis. It is not a fitted event. Specify a
nonzero direction vector; SWCME normalizes it, and the host receives a unit
vector. The direction is expressed in the declared model coordinate frame.
It does not automatically point towards Earth or follow the Parker seed line.
Half width must be in `(0, pi/2]` radians. Axis ratios apply to the unsupported
ellipsoid geometry and do not deform an SSE cap.

For apex distance R and half width lambda, the generating sphere has center
`origin + c*axis`, with `c=R/(1+sin(lambda))`, and radius `a=c*sin(lambda)`.
Along a ray separated by alpha from the axis, the outward surface is

```
R_local = c*cos(alpha) + sqrt(a*a - c*c*sin(alpha)^2)
```

only for alpha <= lambda. The normal is the generating sphere's outward
normal. Its speed is `V_apex*(R_local/R)*(e_r dot n)`. Position, normal and
speed therefore derive from the same self-similar geometry. The tangent-flank
normal speed vanishes. A geometric cap is not automatically a physical fast
shock: SWCME checks the local upstream relative normal speed and fast-mode
speed and creates particle-source patches only where a fast shock exists.

## Application contracts

- Source injection uses `shock_only`/`source`, with the existing DSA population
  normalization and canonical finite source-surface weights. No extra DSA
  acceleration is applied when the mover merely records a geometric crossing.
- Mesh propagation uses schema 4, `background.provider=swcme`,
  `run.intent=shock-propagation`, source disabled, and
  `full_icme`/`resolved_compression`. Combining that resolved compression with
  the DSA source remains rejected to avoid double counting acceleration.
- Outside the cap, canonical SWCME returns ambient Parker/Leblanc primitives.
  No global spherical shock is substituted and no direction is clamped into
  a fabricated flank. Regional widths scale with the local front distance.
- The inner continuation below 1.05 solar radii is ambient. A FULL_ICME epoch
  is accepted only if the entire ejecta and its trailing smoothing transition,
  including the minimum tangent-flank radius `c*cos(lambda)`, stay outside
  that handoff. Wider caps can require a larger launch radius.
- Background generations, owner-cell storage, MPI exchange and DATAFILE
  runtime ownership retain the existing provider-neutral publication path.

## Particle crossing and time steps

`ShockState::MoverGeometry()` transfers geometry, axis, half width and apex
kinematics together during initialization and every native epoch update.
The new header `adapters/shock_geometry.h` defines the AMPS-independent
record. Existing `ExpandingSphericalShock` callers remain source compatible.

Within each substep the apex evolves linearly, as in the preexisting sphere
operator. The SSE generating sphere center translates while its radius grows.
The intersection quadratic uses particle motion relative to that moving
center. Both roots are tested in time order; rear-sphere hits and hits outside
the outward cap are excluded before selecting a crossing. Returned normals
and normal speeds describe the accepted point. Generation de-duplication is
unchanged. Nonlinear DBM kinematics are re-evaluated at the normal runtime
boundaries; reducing the global step tests that temporal approximation.

The pre-step gap to the full generating sphere is a conservative lower bound
on distance to the finite cap. It can select extra substeps outside the cap.
The exact post-step intersection still decides whether a crossing occurred.
A geometric approach is handed to that exact intersection before repeated
gap-halving falls below the minimum numerical substep. Other transport bounds
remain enforced, preventing on-surface source particles from failing merely
because their geometric gap is zero to floating-point precision.

## Mesh coverage and visualization

`examples/sep3d_swcme_sse_mesh_background_20rs_1au.in` retains a Sun-centered
Cartesian box and now defaults to `parker-tube` allocation, with a complete
30-Rs solar neighborhood, a 0.05-AU active tube at 1 AU, and one neighbor-buffer
layer. The 0.03-AU refinement tube is independent of allocation. The input
comments provide exact recipes for switching to `full-domain` (mode plus zero
solar-sphere/tube-radius/buffer sentinels) and back to the corridor defaults.
Full-domain allocation retains both flanks for visualization; corridor allocation
omits distant parts of the CME outside its solar neighborhood and active tube.

The example is a coupling/geometry smoke test. Its minimum cell size is 0.5
solar radii, while the prescribed 0.01-AU-at-1-AU shock transition is only 0.2
solar radii wide at its 20-Rs launch apex, and narrower at the flanks. It is not
resolved for acceleration science. Resolve the shock-normal layer with several
cells, then demonstrate mesh, width and time-step convergence. The declared
memory budget must also accommodate full-domain allocation.

Plot density relative to the unperturbed ambient density, `n/n_ambient(r)`,
or logarithmic density, and overlay the directional front. A linear density
scale dominated by the near-Sun corona obscures compression at 20 Rs.
The propagation CSV's historical `shock_radius_m`/`shock_speed_m_s` columns
mean **apex** distance/speed for SSE. A 1-AU apex crossing does not establish
arrival at an off-axis observer. Direction/width remain in the frozen manifest.

## Build and run

Extract the supplied source overlay at the AMPS root. It contains directly
modified source files, comments and documentation; no patch installer is used.
The overlay also retains the prior runtime DATAFILE ownership and polar-shock
corrections. Refresh any generated application source copies through the usual
AMPS configure/build procedure before linking, then rebuild:

```sh
cd /nobackupp17/vtenishe/Mars1/AMPS
tar -xzf /path/to/SEP3D_finite_SSE_weak_shock_fix_20261002.tar.gz
env MAKEFLAGS="-j16" make -B -C src/models/swcme
env MAKEFLAGS="-j16" make amps
```

Run portable prerequisites from the AMPS root:

```sh
python3 srcSEP3D/test/run_tests.py --group SSE3D --rebuild --output-dir test_output/sse-portable
python3 srcSEP3D/test/run_tests.py --suite standalone --no-build --output-dir test_output/sse-standalone
python3 srcSEP3D/test/test_received_background_mask.py
```

Run the native ten-step coupling check with a fresh output directory:

```sh
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in --test-steps 10 --output-dir test_output/coupled-sep-corona-sse
```

To propagate through the example's configured apex stop, run the normal driver:

```sh
mpiexec -n 4 ./amps --input srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in --output-dir test_output/sse-propagation
```

## Verification and limits

`SSE3D01–09` exercise input validation, 543 canonical geometry comparisons,
mesh fields inside/outside the evolving cap, moving/apex/flank/tangent/rear
particle intersections, real canonical finite source preparation and restart
schema migration, the reported ten-step weak-flank update, and strict
ambient-support/unresolved-layer separation. The native ghost evidence capture now uses AMPS's actual
receive block list and packing mask. Allocated remote centers outside that mask
are not received halo data. `test_received_background_mask.py` compiles the
production mask and evidence-capture functions and verifies all six face
orientations, real stale-data rejection and unmasked receive selection.
This fixture does not claim successful native MPI transport.

`SSE3D05` supplies its injection fixture with a fixed test probe halfway along
the configured finite Parker curve. A generic Earth observer from a different
example can lie outside the active corridor; that is a valid configuration
rejection, not a failed source solver. The fixture uses the actual normalized
curve and coordinate origin, keeps the default `parker-tube` allocation, and
also checks that an opposite-direction observer is rejected. It passes with
`full-domain` allocation as well. See [test/README.md](test/README.md#finite-sse-application-acceptance)
for the setup and [TEST_FIXTURE_FIX_20261002.md](TEST_FIXTURE_FIX_20261002.md)
for the renamed spherical control and runner verification.

New restart payloads use schema 4 and retain finite geometry. Schema-3 files
are read as spheres; ordinary code/configuration identity checks still apply.
Native MPI telemetry also compares shape, axis and half width across ranks.

SSE is prescribed geometry rather than a global MHD simulation. FULL_ICME
still uses phenomenological density/velocity profiles and Parker ejecta B;
there is no flux-rope field or self-consistent CME force balance. The finite
angular edge can truncate those profiles abruptly, so derivatives near that
edge are stencil dependent. This extension establishes finite geometry and
coupling mechanics, not event calibration, divergence-free magnetic ejecta,
resolved shock acceleration or observational agreement.


The October 2 finite-SSE runtime weak-shock correction and its verification
are described in [SSE_WEAK_SHOCK_FIX.md](SSE_WEAK_SHOCK_FIX.md). It resolves
representable weak surface jumps and samples ambient regions without an
irrelevant shock solve; genuine unresolved in-CME states remain explicit.
