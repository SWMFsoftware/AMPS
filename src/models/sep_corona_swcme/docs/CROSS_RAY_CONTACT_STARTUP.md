# Cross-ray startup and field-threaded material contact

Date: 2026-11-10. Status: **bounded geometry and interface diagnostics
implemented/tested; moving-mesh coupling and production selection open**.
This records direct user decisions for the expanded regional closure. It does
not change the legacy `per-ray-lagrangian-piston-v1` fixture or reinterpret its
fixed-HCI contact asset.

## Why the radial startup does not extend unchanged

The retained 3-D ambient corotates in the closed corona. At time zero the
original ellipsoid is stationary in HCI. Its surface is nonspherical and its
centre is displaced from the Sun, so `Omega cross X` is not generally tangent
to that surface. The original supported-ray measurements are:

| Ray | Ambient normal speed (m/s) | Original contact normal speed (m/s) | abs(Bn)/abs(B) |
| --- | ---: | ---: | ---: |
| 82 | -905.1554 | 0 | 0.274351 |
| 83 | +410.1441 | 0 | 0.323732 |
| 106 | -890.4878 | 0 | 0.034639 |
| 107 | +308.6459 | 0 | 0.088526 |

These are actual samples in `test_output/remaining-failures-repair-20261110/contact-topology.csv`.
A zero radial wind does not establish full-vector material impermeability.
An initially ambient plasma, a stationary HCI ellipsoid and a smooth
impermeable boundary are incompatible for this fixture. The original radial
compatibility test remains valid only within its declared reduction.

## Selected angular law and its physical consequences

The user selected corotation at startup, followed by an explicit transition
to the prescribed nonrotating CME motion. The user then selected **T=1200 s**,
matching the existing translational startup duration, and approved the final
angular offset implied by tapering angular velocity.

The reference ellipsoid retains its existing differentiated centre/axis
history. Rotate the WHOLE ellipsoid, including its centre, about the Sun's HCI
z axis. At the owning event start t0 its angle is zero. For
`tau=(t-t0)/T` in [0,1]:

```text
S5(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
omega(t) = Omega_ambient * (1-S5(tau))
theta(t) = Omega_ambient*T*(tau - 2.5*tau^4 + 3*tau^5 - tau^6)
alpha(t) = -Omega_ambient*(30*tau^2-60*tau^3+30*tau^4)/T.
```

After t0+T, omega=alpha=0 and theta=Omega_ambient*T/2. For the retained
`Omega=2.86533e-6 rad/s` and T=1200 s this is 0.001719198 rad, about
0.0985 degrees. No orientation, position or velocity is reset at the join.
The later ellipsoid follows the existing prescribed radial/axis history in
this rotated fixed orientation. This is a new contact identity, not evidence
that an arbitrary observational history should receive the same rotation.

For a fixed reference material label:

```text
X(t) = Rz(theta(t))*X_reference(t)
V(t) = Rz(theta(t))*V_reference(t) + omega_vector cross X(t).
```

`MaterialVelocity` is the velocity of this chosen geometric surface-label
parameterization. Until an ejecta plasma/material closure is supplied, it is
not an additional prescription of the later physical plasma tangential
velocity. Threaded-contact plasma velocity must match on both sides, while
geometry constrains its normal component; tangential relabeling of the surface
cannot change leakage. Actual BG3D-5 plasma matching must distinguish those
contracts. At startup the chosen label motion equals ambient corotation.

The surface normal rotates by the same Rz. Its physical normal speed includes
`(omega_vector cross X).n`. A fixed-HCI ray has a time-dependent direction in
the reference body, so its radial speed/acceleration must differentiate that
rotation as well as centre/axis motion. The implementation uses analytic
implicit derivatives of the same quadratic ellipsoid equation; production
velocity is not a finite difference of large heliocentric positions.

`CorotatingPistonContact::Create(reference,T)` takes T explicitly. It binds the
reference event and T to a new fingerprint, rejects a transition extending
past the declared handoff, and retains the original radial support/status
checks. This low-level interface is implemented; a new runtime model selector,
checksummed expanded input assets, coupled mesh motion and production parser
integration are not yet implemented. The existing input file is untouched.

## Selected magnetic interface

The user selected a **field-threaded ideal material contact**, with compatible
manufactured ejecta data for BG3D-4 boundary tests. Actual two-sided ejecta
matching remains BG3D-5 work. Do not turn the retained nonzero Bn into zero by
projection or claim that this is a closed/draped tangential discontinuity.

For ideal isotropic MHD, with normal n directed ejecta to sheath and physical
surface normal speed wn:

```text
rho_left *(v_left.n - wn) = 0
rho_right*(v_right.n - wn) = 0
v_left = v_right, B_left = B_right, p_left = p_right
T = (p+B^2/(2*mu0))*I - B*B/mu0
[T.n] = 0
[E'_t] = 0, E' = -(v-wn*n) cross B.
```

Density and entropy may jump. The implementation selects the continuous-field
branch everywhere, including its continuous Bn=0 limit. It never switches
branches with a raw floating-point `Bn!=0` check and never divides by Bn. This
does not implement the distinct Bn=0 tangential-discontinuity branch allowing
jumps in tangential field/velocity.

Diagnostics retain signed relative mass flux, full velocity/B/pressure jumps,
normal magnetic jump, Maxwell traction, contact-frame tangential electric
jump and the HCI total-energy interface-flux jump `F_E-wn*E`. They do not
produce an ejecta state, accept a tolerance on the caller's behalf, or close a
force/work budget by editing data. Tangential surface parameterization must
not change physical normal leakage. Actual compatible pressure/field/velocity
and all interface work still need to come from the coupled regional closure.

## Independent tests and exact commands

From `src/models/sep_corona_swcme`, dependency-light shared-library commands:

```sh
make -j16 build/test_cross_ray_contact
build/test_cross_ray_contact
```

These are not native AMPS builds. Four registered cases check:

- `CSWC0628-XRAY-CONTACT-STARTUP`: independent polynomial endpoint/midpoint
  angle/rate/acceleration; actual startup geometry and every original supported
  ray. Maximum startup mass flux is `1.10794e-26 kg/m^2/s`, below the unchanged
  `1e-20` physical bound; max abs(Bn)/abs(B)=0.323732 remains present.
- `...-DERIVATIVES`: independent fourth-order differences of actual material
  positions and fixed-HCI ray radii/speeds through startup, handoff and late
  propagation; tests their distinct derivatives and implicit normal motion.
- `...-INTERFACE`: independent density-jump contact, signed normal leakage,
  shear electric field, oblique rotation, Maxwell stress, pressure work and
  continuous zero-Bn limit.
- `...-REJECTION`: invalid/missing authority, nonpositive or unsupported
  duration, invalid normals/positions/epochs and changed-law fingerprints.

The first manufactured shear expectation mistakenly used 0.05 V/m. Its direct
cross product is `Delta(v_y)*B_x=0.2*0.5=0.1 V/m`. The first FAIL receipt is
retained and the independently derived expected value is corrected; no
production equation or physical gate was changed to pass it. Optimized final
receipts are `contact-rotation-controls-current.json/.log` in the evidence
root. Passing this geometry/diagnostic subset does not close XR1/XR5,
CSWC0628-MASS-TIME, BG3D-4, BG3D-5 or native G1 qualification.
