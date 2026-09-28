# Stage 1: analytic PFSS kernel

The Stage-1 provider implements the current-free coronal field in the shell
`R_sun <= r <= R_b`. It stores photospheric radial-field coefficients in an
explicit real, orthonormal, Condon--Shortley basis. For every non-monopole mode,
the potential coefficients are computed algebraically from the lower radial
field and the condition that the potential is constant at `R_b`. Consequently
`B_theta=B_phi=0` at the source surface by construction, not by clipping.

Associated Legendre functions use a stable degree recurrence. Their first
theta derivative follows the adjacent-degree identity and the second follows
the associated-Legendre differential equation. Radial factors and their first
two derivatives are analytic. `SphericalField::derivative` therefore exposes
all nine partial derivatives without finite differencing the physical field.
Cartesian queries transform the final vector only; the harmonic authority stays
spherical.

The monopole is never part of a PFSS solve. `Reject` fails when it exceeds the
declared tolerance, while `Remove` records and removes it during map projection.
The heat-kernel operation applies
`exp[-l(l+1)/(l_a(l_a+1))]` to a returned copy. `Scaled` is also pure, so the
unscaled coefficient authority cannot be mutated or scaled twice.

Field lines use fourth-order Runge--Kutta integration of `dx/ds=+/-B/|B|` in
Cartesian coordinates. Classification is based on reaching a physical radial
boundary, never on an AMR cell. PFSS-to-`R_b` and future composite-to-`R_i`
classifications remain separate in `TopologyPair`; only the composite result
routes open-wind versus closed-hydrostatic plasma ownership.

Map projection requires explicit solid-angle weights. Tests cover a dipole,
an individual non-axisymmetric mode, zero signed flux, source-surface radiality,
step-converged topology, map round trip, heat-kernel transfer, deterministic
direct evaluation, immutable scaling, and independent topology routing.
