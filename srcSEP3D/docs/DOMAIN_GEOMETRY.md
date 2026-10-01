# Domain geometry, solar neighbourhood and finite field-line corridor

The complete example is `examples/sep3d_analytic_parker_corner_sphere.in`.
These controls extend the existing input grammar additively. Older decks keep
their Sun-centered cube, source-shell radial anchor and corridor-only mask.
The Sun remains at the heliocentric coordinate origin; the Cartesian root
box moves around it. No field, CME, observer or coordinate frame is translated.
The supplied example places the Sun near an x-y corner and on the z midplane;
its complete near-Sun sphere and equatorial corridor extend to both sides of
z=0. The earlier three-axis corner layout remains an explicit alternative.

## Input controls

```ini
[domain]
inner_radius_m = 6.957e8
outer_radius_mode = field-line-endpoint
box_geometry = field-line-xy-corner-cube
corner_direction_x = 0
corner_direction_y = 0
corner_direction_z = 0
corner_margin_m = 1.495978707e9

[parker_spiral]
start_mode = explicit
initial_x_m = 6.957e8
initial_y_m = 0
initial_z_m = 0
end_radius_m = 1.6455765777e11
length_m = 0

[mesh.solar]
anchor = photosphere
enabled = true
surface_cell_size_m = 3.4785e8
transition_outer_radius_m = 1.3914e10
profile = smoothstep
exponent = 1

[mesh.active_region]
mode = parker-tube
solar_sphere_radius_m = 2.0871e10
```

This is a fragment, not a replacement for a complete deck. Existing corridor
width, refinement, halo, provider, species, source and output controls remain
explicit in the supplied example. All new corner fields must be present when
that layout is selected. `outer_radius_m` and `preset` are inactive in endpoint
mode; the resolved endpoint distance enters the immutable configuration.

`end_radius_m` is a heliocentric distance in metres. It owns the physical outer
cutoff and finite Parker line length, calculated with the canonical arc-length
law. `length_m=0` requests derivation; a positive matching normalized length is
accepted for typed/parser round trips, while a conflicting length is rejected.
`point_count` controls diagnostic sampling and cannot change the root bounds.

For x/y, each `corner_direction` is -1, +1, or 0. Positive means that axis extends into
the box from a lower corner face; negative selects an upper face. Zero selects
the sign from the endpoint, with a whole-curve fallback for a zero component.
The Sun is inset by at least `solar_sphere_radius_m + corner_margin_m` from
the chosen x/y faces. The complete sphere fits: an exactly corner-centered Sun
would permit only part of the sphere and is deliberately not this layout.

| `box_geometry` | Position of the Sun relative to the root box |
|---|---|
| `sun-centered-cube` | Centered in x, y and z; original input default. |
| `field-line-corner-cube` | Near a selected corner in all three axes; z uses its own direction selection. |
| `field-line-xy-corner-cube` | Near a selected x-y corner, centered between the two z faces. |

In the x-y mode, `corner_direction_z=0` is required and denotes a centered z
extent. A nonzero value is rejected rather than silently ignored. The two z
faces are `origin_z_m - side/2` and `origin_z_m + side/2`. The production
heliocentric contract keeps the Sun at `(0,0,0)`, so these bounds are symmetric
around zero. Centering the Sun between the faces is a change to the box bounds,
not a translation of the Sun, the Parker field or any observer.

The common cube side begins at `end_radius_m + 2*(sphere_radius + margin)`.
Whole-curve/cross-section bounds can increase the side or an axis inset when
a winding bends towards another face. The enclosure uses the analytic bound
on `|dx/dr|` over radial intervals, independent of plotted points. Long,
tightly wound curves can need more space than their endpoint suggests.
The margin is a physical box clearance, not a halo-block count. Choose it to
support the intended boundary stencils and inspect the actual AMR output.
For x-y mode, the side also satisfies
`side >= 2*(max(sphere_radius, abs(curve_z_min), abs(curve_z_max)) + margin)`;
the curve bounds already include the full corridor width. A tilted or polar
selected line therefore enlarges the symmetric z extent, and the common cube
side, instead of clipping its cross-section or moving the Sun off the midplane.

For the complete supplied 1.1-AU example, the resolved bounds in metres are:

```text
Sun center: (0, 0, 0)
minimum: (-2.2366978707e10, -1.86924636477e11, -1.04645807592e11)
maximum: ( 1.86924636477e11,  2.2366978707e10,  1.04645807592e11)
```

Thus z spans approximately `[-0.699514, +0.699514] AU`; the Sun is at the exact
z midplane and near the negative-x/positive-y corner. The solar allocation
sphere is 30 R_sun, with the configured extra clearance at those x/y faces.

## Active allocation and physical boundaries

The physical core is the union of the finite corridor and the complete solar
neighbourhood. Exact sphere/AABB intersection retains every intersecting leaf,
including face tangencies. Existing face/edge/corner halo layers and cavity
checks then apply to that union. AMPS independently removes leaves wholly
inside the fixed 1-R_sun photosphere and preserves fractional cut-cell measures.
The allocation sphere is not another absorbing boundary or a new solar radius.

The sphere must reach the configured Parker/CME source radius and lie strictly
inside the endpoint radius. This keeps the spherical and corridor components
connected. `solar_sphere_radius_m=0` disables the spherical addition in a
legacy centered corridor; a corner layout requires it. This geometry change
does not change `inner_radius_m`, the transport/source cutoff, or the selected
background provider and its physical validity. The supplied corner/sphere
example has a 1-R_sun absorbing inner boundary/Parker source and a 1.1-AU
outer transport cutoff and line endpoint. Its 30-R_sun allocation sphere and
20-R_sun refinement transition are separate geometry controls. The SWCME
event still launches at 20 R_sun: `start_mode=explicit` permits the line to
start at the photosphere without changing the CME launch event. The canonical
SWCME Parker source, radial field cross-check and on-line observer coordinates
are updated consistently with the photospheric start. One-AU reference radii
still normalize fields, turbulence and tube widths; they are not boundaries.

## Radial coarsening

Let `r0` be R_sun for `anchor=photosphere`, or `inner_radius_m` for the legacy
`anchor=source-shell`. Let `r1=transition_outer_radius_m`, and clamp
`x=(r-r0)/(r1-r0)` to [0,1]. The radial target is

`h(r)=surface_cell_size_m + f(x)*(global_cell_size_m-surface_cell_size_m)`.

The profiles are `f=x` for linear, `f=x^exponent` for power-law, and
`f=[x*x*(3-2*x)]^exponent` for smoothstep. For the latter two profiles, a larger
positive exponent retains fine resolution farther out; a smaller exponent
coarsens sooner. Linear ignores the exponent. At/below the anchor, the surface
target is exact; at/beyond the transition, the radial target is global.
Tube refinement may request a finer target there: the combined law takes the
minimum of radial, tube and global requests, bounded by the mesh minimum.

The physical photosphere uses this same AMPS `localResolution` callback.
The portable octree also tests the nearest box point to the Sun so an off-grid
solar origin cannot hide a small radial refinement region. Maximum-level
feasibility and dry-run volume sampling use the actual shifted cube side and
coordinates, rather than the legacy `2*outer_radius_m` estimate.

## Commands and verification

From the AMPS root, inspect without allocating the mesh:

```sh
./amps --input srcSEP3D/examples/sep3d_analytic_parker_corner_sphere.in --dry-run
```

Run the existing native initialization suite with the new deck:

```sh
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_analytic_parker_corner_sphere.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/corner-sphere/native.json --artifact-directory test_output/corner-sphere/artifacts/
```

Native SCCM3D04 checks installation and allocation of the configured active
plan; it remains one generic host case. The application `test/run_tests.py
--all` additionally selects DOM3D01–04 and CFG3D12. DOM3D04 covers all x/y
corners, automatic selection, tilted/polar corridor enclosure and reflected-z
mesh/active-plan symmetry. CFG3D12 parses and fingerprints both corner modes.
The shared-model/native
aggregate driver does not run the application-only portable geometry group.

Dry-run and native startup report geometry, minimum/maximum box coordinates,
endpoint distance, sphere radius and the radial anchor. Geometry controls and
resolved bounds enter the physics fingerprint, so incompatible restarts fail.
Legacy decks retain their existing fingerprint when these controls are unused.
Portable geometry/AMR and runner tests are executable without a full checkout;
real parser/provider and multi-rank AMPS verification need the enclosing legacy
shared-model sources and configured AMPS build.
