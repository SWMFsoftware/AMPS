# Material-boundary Maxwell–Jüttner flux design

Status: the **distribution and material-grid choice are adopted by the user**.
The sampler/normalization/production bridge described here is a derived design,
not an implemented or qualified source. Track implementation in
[STATUS.md](../../STATUS.md). No classical tail truncation is authorized.

## Required physical inputs and measure

Specify proper species number density `n0` [m^-3], proper temperature `T` [K],
rest mass `m` [kg], the full bulk velocity `U` [m/s], and an oriented instantaneous
boundary area/normal in the simulation frame. The material boundary has
`v_grid=U` at the source event. The density is not silently taken from an
electron-density field or converted to proper density without a provider contract.

The isotropic rest distribution is the standard Maxwell–Jüttner distribution
with probability measure `f0(p') d^3p'`, proportional to
`exp[-gamma'/theta] d^3p'`, where `theta=kB*T/(m*c^2)`. This is a volume
distribution before applying the boundary-crossing weight. The proper-frame
definition and drift transformation must stay explicit; see
[Swisdak (2013)](https://arxiv.org/html/1305.5226v1) and
[Zenitani (2015), equations 22–28](https://arxiv.org/html/1504.03910v1).
In particular, `d^3p/E` is invariant, so a simple boost of a volume sample
omits the momentum-space Jacobian in a simulation-frame density/flux measure.

This derivation describes **full-particle** crossings of a moving plane.
The srcSEP transport state is gyroaveraged and retains only parallel and
perpendicular momentum magnitudes. A tilted material end cap can receive a
gyrophase-dependent crossing contribution that this reduced state cannot
represent. The native bridge therefore needs a derived guiding-center boundary
matching law, or an explicit domain where the transformed cap normal is
parallel to the proper-frame magnetic field and the reduction is valid. Merely
projecting a general full-flux sample onto two momenta does not qualify that
boundary condition. Proper-frame B and simulation-frame B must also remain
distinguished when applying a Lorentz boost.

## Derived material-boundary sampling law

The following algebra is a design derivation from that measure, rather than a
claim that a volume-loading algorithm from a reference is already a flux sampler.
Let `n` be the unit normal pointing into the simulated domain, `Gamma` the bulk
Lorentz factor, `Uhat=U/|U|`, and use the continuous `U=0` limit when appropriate.
For a boundary normal velocity `vgn=v_grid dot n`, define

```text
C = Gamma*(U dot n - vgn),
a = n + (Gamma-1)*(Uhat dot n)*Uhat - Gamma*vgn*U/c^2.
```

Multiplying relative normal speed by the momentum Jacobian gives exactly

```text
(E/E') * max[(v-v_grid) dot n,0]
    = max[C + (c^2/E') * a dot p',0].
```

For the selected material boundary `v_grid=U`, `C=0` and
`|a|=sqrt[1-(U dot n/c)^2]`. In terms of `x=gamma'-1`, the rest radial density
before flux weighting is

```text
g(x) dx proportional to (1+x)*sqrt[x*(x+2)]*exp(-x/theta) dx.
```

The angular integral of the positive crossing factor adds
`c*|a|*sqrt[x*(x+2)]/[4*(1+x)]`. Therefore the normalized *flux* radial law is

```text
q(x) proportional to x*(x+2)*exp(-x/theta), x>=0.
```

This is an exact mixture: choose `Gamma(shape=2,scale=theta)` with probability
`1/(1+theta)`, or `Gamma(shape=3,scale=theta)` with probability
`theta/(1+theta)`. Each integer-shape gamma can be sampled with respectively two
or three independent open-unit exponential variates. No momentum cutoff,
low-temperature switch to another law, or radial rejection approximation is
needed. Evaluate mixture probabilities without overflowing `1+theta`.

About the axis `a/|a|`, sample rest-frame cosine `mu=sqrt(U_open)` and a uniform
azimuth. Construct the full three-component rest momentum with magnitude
`m*c*sqrt[x*(x+2)]`, then perform the full Lorentz boost if a simulation-frame
momentum is needed. Keep the original rest momentum for any explicitly
plasma-frame native storage; avoid a numerically destructive boost/inverse-boost
round trip. Boundary acceptance, magnetic projection and source normalization
must use their declared frames.

This simplification is **material-boundary only**. A nonmaterial boundary has
`C!=0` and a different angular/radial selection law. It must not call this kernel
as an undocumented approximation.

## Physical normalization and numerical obligations

The number flux per instantaneous simulation-frame area/time is

```text
F = n0*|a|*<v'>/4,
<v'> = 2*c*theta*(1+theta) / [exp(1/theta)*K2(1/theta)].
```

`K2` is the modified Bessel function. Its exponentially scaled form, or an
independently bounded positive integral, is required; separately forming
`exp(1/theta)` and `K2` fails at low temperature. A source implementation must
record its numerical tolerance, quadrature/tail bound or special-function
reference and failure status. An unrepresentable result is explicit
`Unresolved`, not zero flux or an alternate sampler.

As `theta` tends to zero,
`F/(n0*|a|)` approaches `sqrt[kB*T/(2*pi*m)]`. This is the classical Maxwellian
flux limit; it does not authorize using a classical distribution at finite
temperature. For `U=0`, the limiting inward normal velocity is flux weighted,
with mean `sqrt[pi*kB*T/(2*m)]`, normal second moment `2*kB*T/m`, and each
transverse second moment `kB*T/m`.

The source count additionally requires physical face area, proper density,
statistical particle weight and time integration of the flux. A uniformly
distributed birth phase is exact only for a declared constant source rate
within its interval. Time-varying rate/area/temperature requires its explicit
integration/sampling policy and convergence evidence. The sampler alone does
not settle that native bridge.

## Independent qualification gates

- Verify the invariant-measure/boosted moving-plane flux against direct
  multidimensional quadrature, including normals parallel and transverse to U.
- Check the exact two/three-gamma mixture CDF and moments with predeclared
  ensemble sizes, seeds and statistical limits; retain separate radial and
  angular controls.
- Verify low-temperature flux and moments without swapping the implemented
  law; test finite relativistic temperatures and bulk velocities separately.
- Check subluminal velocities, material-boundary orientation, finite positive
  inputs, numerical representability and atomic failure before RNG/source
  counters or particle allocation are committed.
- Qualify actual keyed count/phase/momentum sampling, native rest-frame
  storage, material face geometry, remaining birth duration, MPI ownership,
  normalization and checkpoint continuation. No source-only sampler test
  establishes those production seams.
