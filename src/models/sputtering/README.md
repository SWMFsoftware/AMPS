# AMPS surface sputtering model

This directory contains the ice sputtering-yield and emitted-speed models used
by AMPS, plus a dependency-light emission and conservation-accounting layer.
The new source layer does not allocate particles, select an AMR owner, move a
particle, or reduce MPI sampling buffers. Those remain application/PIC tasks.

## Implemented physical interfaces

`Sputtering::GetSputteringYield()` returns the dimensionless number of emitted
products per incident projectile for the configured ice model. The retained
models are:

- `Ice::Rubin`: H2O production for H+, O ions, and O2+;
- `Ice::Teolis`: O2 and H2O production for the same projectile families.

Both tables use projectile kinetic energy

\[
E_{\rm eV}=\frac{m v^2}{2e}
\]

with mass in kg and speed in m/s. Positive tabulated yields are interpolated
linearly in `log10(E)`–`log10(Y)`. Out-of-domain and unsupported projectile
queries return zero; unsupported product identities use the established AMPS
fatal-configuration path. The table comments identify Martin Rubin and
Benjamin Teolis as the model sources, but this repository does not retain an
exact publication, version, digitization record, or uncertainty covariance for
either table. Consequently the tests qualify implementation and interpolation,
not independent laboratory provenance.

`Ice::GetSputteringSpeed()` samples uniformly in `ln(v)` over the retained
finite ranges and returns a relative importance weight. Its H2O and O2 energy
laws and binding-energy parameters are documented beside the implementation.
No evidence file in this directory identifies their original publication, so
those parameters must not be presented as independently revalidated data.

## Emission-angle selector

[`Sputtering_Source.h`](Sputtering_Source.h) provides two runtime-selectable
outward-hemisphere laws:

| Input selector | Probability density | Inverse CDF for `mu=cos(theta)` |
|---|---:|---:|
| `cosine-knudsen` | `p(Omega)=mu/pi`, `p(mu)=2mu` | `mu=sqrt(u)` |
| `uniform-hemisphere` | `p(Omega)=1/(2pi)`, `p(mu)=1` | `mu=u` |

`cosine-knudsen` is the default. The cosine factor is the projected-flux
factor for an isotropic velocity population crossing a plane; it is not an
ice-specific empirical fit. `uniform-hemisphere` is a parameter-free control
law already supported by the common AMPS surface-source layer. Neither law
introduces an energy-angle correlation.

Example:

```cpp
using Sputtering::Source::EmissionAngleLaw;

const EmissionAngleLaw law =
    Sputtering::Source::ParseEmissionAngleLaw("cosine-knudsen");
const std::array<double, 3> direction =
    Sputtering::Source::DirectionFromUnitSamples(
        outward_normal, random_u, random_phi, law);
```

The normal and returned direction use the caller's Cartesian frame. The normal
is normalized internally, and the returned vector points outward. Unknown
selector names, invalid enum values, invalid variates, and zero/nonfinite
normals produce typed `std::invalid_argument` failures. Future laws should add
one enum value, its stable selector name, an explicit normalized PDF/inverse
CDF, analytic moments, and independent tests. They must not silently reuse a
different law.

## Source and statistical-event accounting

For incident number flux `Gamma` in projectiles m^-2 s^-1, physical area `A`
in m^2, and product yield `Y` in products/projectile,

\[
R_{\rm incident}=\Gamma A,\qquad
R_{\rm product}=\Gamma A Y.
\]

If a numerical emission event represents `W` physical products, its mean event
rate is `lambda=R_product/W`; over physical interval `dt`, the expected event
count is `lambda*dt`. `ComputeSourceRateAccounting()` evaluates these equations
with finite-domain and overflow checks. `EffectivePhysicalWeight()` evaluates
`W=W0*c` from base statistical weight `W0` and individual correction `c`.

For variable-weight events, a nominal or arbitrary correction is not a valid
normalization. The caller must either use a qualified expectation for event
scheduling or close realized event weights with `ProductLedger`.

## Product conservation ledger

`ProductLedger` is bound to one nonnegative product-species identity. For an
event representing `W` physical products of mass `m` and velocity `v`, it
accumulates

\[
N\mathrel{+}=W,\quad M\mathrel{+}=mW,\quad
\mathbf{P}\mathrel{+}=mW\mathbf{v},\quad
K\mathrel{+}=\tfrac12mW|\mathbf{v}|^2.
\]

Dividing by elapsed physical time produces number rate in s^-1, mass rate in
kg/s, momentum-transfer rate in newtons, and kinetic power in watts. Velocity
must be expressed in one caller-declared inertial or body-fixed Cartesian
frame. Momentum is positive along the supplied velocity; for the angular API,
the surface normal establishes the positive outward direction.

The ledger contains emitted-product translational mechanics only. It does not
claim closure of projectile recoil, surface binding/thermal energy, lattice
momentum, electromagnetic work, chemistry, or internal degrees of freedom.
Species mismatches and invalid values are typed errors; updates are
transactional so a rejected event does not partially alter a ledger.

## Current integration boundary

`Sputtering.h` exposes the source layer, so every application already compiling
the sputtering module can opt into it without another library. Existing
application injection behavior is unchanged. In particular, legacy
application-local angle generators have not been silently replaced because
doing so would change their physical distribution. Migration requires an
application-specific comparison and explicit selection of one documented law.

PIC particle allocation, AMR ownership, time-of-event sampling, mover
insertion, application sampling buffers, and MPI reduction are not qualified
by the dependency-light suite. A successful source-layer test must not be
reported as evidence for those native integration paths.

## Tests

The focused runner is:

```text
python3 src/test/tests/models/sputtering/run_test.py --context <context.json>
```

The global AMPS test runner supplies that resolved context. The suite compiles
the production yield, speed, angle, rate, and ledger code with a narrow PIC
contract stub and checks:

- Rubin/Teolis endpoints, thresholds, product identities, and interpolation;
- H2O/O2 speed-law importance weights and finite-window normalization;
- selector/default/typed-failure behavior;
- analytic cosine and uniform-hemisphere angular moments;
- `flux * area * yield` physical rates and statistical event expectations;
- product number, mass, momentum, kinetic-energy, and time-normalized ledgers;
- joint planned-rate and sampled-ledger closure for both angular laws.

These tests are deterministic and use independent equations or pinned retained
table points. They do not generate expected values by calling a second copy of
the production algorithm.

The registered portable user-facing commands are:

```text
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-SPUTTERING-001 --order normal
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-SPUTTERING-001 --order reverse
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-SPUTTERING-001 \
  --order shuffle --shuffle-seed 707
```

The suite is a direct dependency-light compile, not a native AMPS rebuild. A
future native integration test must use the repository's required clean build,
configuration, and top-level `-j16` workflow; a source-only compile must not be
reported as native evidence.
