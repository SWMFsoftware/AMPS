# Phase P: Three-Dimensional Transport Cores

Phase P implements the AMPS-independent Parker and focused-transport steps.
The production AMPS adapter translates particle and cell records into these
types; it does not contain a second transport equation. All coordinates and
coefficients use SI units.

## Tensor Parker step

For parallel-only diffusion, the spatial diffusion tensor is

\[
\boldsymbol{\kappa}=\kappa_\parallel\mathbf b\mathbf b,
\]

where \(\mathbf b=\mathbf B/|\mathbf B|\). The implemented Itô SDE is

\[
d\mathbf X=[\mathbf U+\nabla\!\cdot\boldsymbol{\kappa}]dt
 +\sqrt{2\kappa_\parallel}\,\mathbf b\,dW,
\qquad
dp=-\frac{p}{3}(\nabla\!\cdot\mathbf U)dt.
\]

The complete tensor divergence is evaluated as

\[
\nabla\!\cdot(\kappa_\parallel\mathbf b\mathbf b)=
\mathbf b\frac{d\kappa_\parallel}{ds}+
\kappa_\parallel[(\mathbf b\!\cdot\nabla)\mathbf b+
\mathbf b(\nabla\!\cdot\mathbf b)].
\]

Thus `ParkerLocalState` requires the field-aligned coefficient derivative,
field-line curvature, and `divBhat`. Omitting either geometry term would change
the Fokker–Planck operator and destroy the uniform-density equilibrium for a
nonuniform tensor. The frozen-`divU` momentum characteristic is integrated
exponentially, preserving positive momentum. `COEF3D03–05` verify tensor rank,
second-order agreement with a numerical tensor divergence, and rejection of
invalid coefficients. `PRK3D01–08` cover moments, advection, rotations,
nonuniform equilibrium, cooling, first passage, a PDE Green-function baseline,
and every named timestep limiter.

## Focused transport step

The gyrotropic state is `(positionM, momentumKgMPerS, mu)`. For one immutable
local background sample, the deterministic coefficients are

\[
\dot\mu=\frac{1-\mu^2}{2}
\left[v\nabla\!\cdot\mathbf b+
\mu(\nabla\!\cdot\mathbf U-3\mathbf{bb}:\nabla\mathbf U)\right],
\]

\[
\frac{d\ln p}{dt}=-\frac12\left[(1-\mu^2)\nabla\!\cdot\mathbf U+
(3\mu^2-1)\mathbf{bb}:\nabla\mathbf U\right],
\]

and spatial streaming is \(\dot{\mathbf X}=\mathbf U+\mu v\mathbf b\).
Pitch-angle scattering solves the conservative operator
\(\partial_\mu(D_{\mu\mu}\partial_\mu f)\), so the Itô pitch SDE contains
the required drift \(\partial D_{\mu\mu}/\partial\mu\) and noise
\(\sqrt{2D_{\mu\mu}}dW\).

The algorithm is a symmetric deterministic–stochastic–deterministic split:

1. Advance focusing/flow pitch drift and momentum for half a step.
2. Apply one keyed stochastic pitch kick. The default scalar Milstein
   correction is
   \(\tfrac12D'[(dW)^2-dt]\); reflecting Euler–Maruyama remains an explicitly
   fingerprinted comparison option.
3. Mirror-fold any finite overshoot at \(\mu=\pm1\), implementing the same
   zero-flux boundary as the Fokker–Planck equation without endpoint clipping.
4. Advance the second deterministic half-step and stream with midpoint pitch
   angle and speed.

The pitch-averaged momentum equation reduces to Parker cooling. `FTE3D01–07`
verify ballistic motion, focusing and mirroring, Legendre eigenmode rates,
bounded strong scattering, momentum characteristics, the strong-scattering
Parker limit, and bitwise identity of the zero-perpendicular hooks.

## Named timestep limits

`SelectTimeStep` reports the exact minimum of requested, cell-crossing,
diffusion, focusing, cooling, fractional-field-variation, shock-crossing, and
snapshot-boundary limits. The fixed candidate order is part of the diagnostic
schema. A value below `minimumTransportSubstepS` returns `StepUnderflow` with
the actual value and limiter; it is never replaced by the minimum.

## Reproducible random streams

`KeyedRandomStream` hashes campaign seed, stable particle ID, step, substep,
semantic purpose, and draw counter. It contains no shared engine and no rank or
thread identifier. One normal draw always consumes two counters. Consequently:

- particle-list order and MPI/OpenMP ownership cannot change a history;
- restart serializes one small counter state per active purpose;
- adding a future random purpose cannot shift Parker, pitch, or source draws.

`RNG3D01–03` exercise worker partitioning, reversed particle order, and an
otherwise unused future-physics stream bit-for-bit.

## Deliberately unavailable physics

Both cores require perpendicular diffusion and explicit drift velocity to be
exactly zero. A nonzero value returns `InvalidInput`. This is a release guard:
those terms require their own coefficient, timestep, boundary, and validation
campaign and are not represented by silent no-op branches.
