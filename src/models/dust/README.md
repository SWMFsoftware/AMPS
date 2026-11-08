# AMPS dust model

This directory contains the PIC-coupled electrically charged dust model and
two dependency-light physical contracts:

- `DustCore.h` implements logarithmic size mapping, spherical grain mass,
  isolated-sphere potential, Lorentz acceleration, isotropic-direction
  normalization, power-law size normalization, and weighted sampling moments.
- `DustSplitMerge.h` implements deterministic dust-specific statistical
  splitting and merging without changing the generic PIC split/merge code.

The dependency-light routines are production APIs used by `Dust.h` and
`Dust.cpp`; they are not independent replacement physics models.

## Units and grain state

Radius is in metres, physical grain mass in kilograms, grain charge in
coulombs, velocity in metres per second, electric field in volts per metre,
magnetic field in tesla, and acceleration in metres per second squared.  A
dust macroparticle represents its species-global base weight multiplied by its
positive `individual_stat_weight_correction`.  Dust split/merge never falls
back to unit weight when the individual-weight build capability is disabled.

The isolated spherical potential is

```text
phi = q / (4 pi epsilon_0 a),
```

and the point-grain Lorentz acceleration is

```text
dv/dt = (q/m) (E + v cross B).
```

Neither formula includes plasma shielding, drag, gravity, or radiation
pressure.  Those effects require their own enabled closures.

## Size groups

`SplitMerge::SizeGroups(a_min, a_max, count)` defines equal-width intervals in
`ln(a)`.  Group `i` is

```text
[a_min exp(i Delta), a_min exp((i+1) Delta)),
Delta = ln(a_max/a_min)/count,
```

with `a_max` owned by the final group.  Values outside `[a_min,a_max]` are
invalid split/merge inputs rather than clamped candidates.

Two particles may be merged only when all of these conditions hold:

1. individual statistical weights are enabled, finite, and positive;
2. particles have the same native ownership-cell key and compiled species;
3. radii map to the same declared size group; and
4. radius, per-grain mass, and per-grain charge are exactly identical.

The same-group rule is necessary but not sufficient.  Exact intrinsic equality
is intentional because continuous radius changes charging, drag, radiation,
and sampled-size physics.  A tolerance-based merge could silently change those
closures and make mass/charge conservation tolerance-dependent.

## Deterministic split and merge

`Split` copies one complete grain state into `n >= 2` representatives and sets
each individual correction to `w/n`.  It exactly preserves represented number,
mass, charge, momentum, kinetic energy, size group, and all current sampling
moments up to floating-point summation.

`Merge` requires at least three admissible particles and produces two.  If
`W=sum(w)` and `u=sum(w v)/W`, both outputs receive `W/2`; their velocities are
`u +/- d e_j`.  Here `d^2=sum(w |v-u|^2)/W`, and `j` is the velocity component
with the largest variance, with x-before-y-before-z tie breaking.  This
deterministically preserves represented number, mass, charge, vector momentum,
and scalar kinetic energy.  It also preserves weighted velocity and electric
current.  It does not generally preserve the pressure tensor, higher velocity
moments, or `sum(w |v|)`; callers requiring those quantities must retain the
original population.  Passing `require_sampled_speed_conservation=true`
returns the typed `kSampledSpeedUnrepresentable` outcome and no outputs when
the conservative two-particle state cannot also close that raw sampled moment.
All validation precedes construction of output states, so a rejection is
transactional.

Typed outcomes distinguish invalid configuration/particles, disabled
individual weights, different owner/species/group/intrinsic grain, and an
insufficient population.  A future PIC list adapter must validate first, then
copy all particle extension bytes and update list links in the same native AMR
cell.  Generic `PIC::ParticleSplitting` remains unchanged because it does not
reconstruct dust radius, mass, or charge and uses a randomized merge direction.

## Spherical injection

`DustInjection__Sphere` currently draws an isotropic surface normal from
independent Gaussian components, assigns radial speed, consumes a residual
mass reservoir, applies AMR ownership and local time-step acceptance, creates a
particle, samples surface mass flux, and dispatches the configured mover.
The source state uses `x=x0+R*n`, so translated internal spheres retain their
configured centre.  Temporary particle bytes are initialized, a missing AMR
owner is rejected before dereference, and `Injection::ResetMassReservoir()` /
`GetMassReservoir()` provide an explicit per-rank test and initialization
boundary for the signed residual-mass ledger.

The algebraic direction, translated-sphere state, and grain-mass laws are
tested by `CORE-DUST-001`.  The separate
`CORE-DUST-INJECTION-NATIVE-001` profile compiles DUST with individual weights,
constructs a one-rank 3-D AMR mesh and translated spherical boundary, invokes
the production injector/mover twice, and observes the native owner, particle
record, sampled mass rate, and signed residual ledger.  The deterministic
source requests half a represented grain per interval: the first interval
emits one grain and leaves a negative half-grain residual, and the second emits
none and returns the reservoir exactly to zero.  This qualifies bounded
one-rank insertion; multi-rank source arbitration, migration, restart, and MPI
reduction remain outside T7.

## Focused tests

Run the registered dependency-light suite from the repository root:

```bash
python3 src/test/run_tests.py --test CORE-DUST-001 --jobs 1 \
  --output-dir /tmp/amps-core-dust
```

The suite uses literal or independently derived analytic fixtures.  A native
AMPS rebuild is not part of this command.  Run the configured integration case
with:

```bash
python3 src/test/run_tests.py --test CORE-DUST-INJECTION-NATIVE-001 --jobs 1 \
  --output-dir /tmp/amps-core-dust-injection
```

That exclusive runner performs the repository `AGENTS.md` native preflight,
clean Config.pl selection, and top-level `-j16` build before the one-rank run.

