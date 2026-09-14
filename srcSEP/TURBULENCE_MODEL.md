# Alfvén-turbulence ownership and evolution contract

Steps 11 and 12 retain self-consistent Alfvén turbulence as a physics
subsystem; it is not a fourth particle mover. The exactly three production
movers remain `parker`, `fte-dmumu`, and `fte-mfp`.

## Source and authoritative state

`util/sep_turbulence_core.*` defines five explicit source modes:

| Source | Mutable here? | Authoritative representation |
|---|---:|---|
| `prescribed` | No | selected prescribed state |
| `self-consistent-integrated` | Yes | segment-integrated `E+`, `E-` [J] |
| `self-consistent-spectral` | Yes | branch-major `E±(k)` [J per log bin] |
| `swmf-read-only` | No | immutable imported state |
| `swmf-initial-then-local` | After one handoff | selected local representation |

Integrated and spectral arrays cannot be simultaneous mutable owners. In
spectral mode, branch-integrated values are derived by summing bins. In
integrated mode, spectral storage must be absent. An SWMF-to-local handoff is a
one-time operation and records its epoch and source checksum.

Wave density is derived as `W±=E±/V` [J/m³], amplitude as
`deltaB²=2*mu0*(W+ + W-)`, and cross helicity as
`sigma_c=(E+ - E-)/(E+ + E-)`. Coefficient views derive resonance, branch,
`Dmumu`, parallel mean free path, and parallel spatial diffusion from this same
state; they do not own a second copy.

## Driver order and stability

`SEP::Turbulence::Advance` is the context-neutral driver for standalone and
coupled state adapters. Its fixed order is:

1. shock and other external wave sources;
2. particle-wave energy exchange;
3. finite-volume advection with explicit CFL subcycling;
4. branch reflection;
5. cascade and physical dissipation;
6. authoritative/derived representation synchronization; and
7. diagnostics, ledger accumulation, and checkpoint phase reset.

The standalone loop installs the same active source, operator, and coefficient
configuration. Prescribed and `swmf-read-only` sources bypass every mutating
operator but remain available for synchronization and diagnostics. An invalid
or non-positive production CFL limit is fatal; it is never replaced by the
remaining particle timestep.

Reflection is an internal signed branch transfer and preserves total energy.
Particle and shock withdrawals are limited by available branch energy, with
the requested/applied difference reported as a limiter correction. The
integrated reference cascade uses an implicit positive turnover sink.

## Boundary and conservation ledger

Each field-line end owns one restartable boundary object:
`specified-incoming-energy`, `specified-incoming-flux`,
`transparent-outflow`, or `fixed-reservoir`. Energy values use joules and flux
values use inward watts. Non-transparent negative or non-finite values are
rejected.

Every update reports signed initial/final energy, inner and outer boundary
exchange, shock source, reflection transfer, cascade transfer, physical
dissipation, particle exchange, limiter correction, remap correction, and
closure residual. Reflection and internal cascade transfer cancel from the
total. Material positivity corrections increment the limiter diagnostic.

## Remap and restart

`RemapConservatively` overlaps old and new cells in physical arc length. It
conserves each integrated branch and, in spectral mode, every individual bin on
the covered domain. It increments the field-line generation and reports any
domain or roundoff correction.

The deterministic versioned checkpoint contains source, representation,
operators, boundary objects, spectral grid, phase, epoch, field-line
generation, campaign seed, completed step, pending shock/particle exchange,
handoff epoch/checksum, provenance, accumulated ledger, geometry, and all
authoritative energies. Deserialization validates ownership and dimensions.

## Deterministic concurrency policy

`util/sep_reproducible_reduction.*` supplies worker-local append-only buffers.
Contributions are keyed by field line, segment, branch, particle, timestep, and
purpose—never MPI rank or thread. After workers join, the reducer sorts by the
complete key and accumulates with long-double intermediates. MPI rank blocks
must be gathered and passed through the same canonical reducer. `MPI_Allreduce`
is reserved for genuinely additive accumulators; non-additive authoritative
state uses owner/gather synchronization.

Warnings and limiter counts are integer atomics. Random streams are keyed by
campaign seed, particle ID, timestep, and purpose so optional diagnostics cannot
advance a mover stream. `EvidenceHashHex` hashes the canonical reduced result;
PAR01–PAR05 require the same hash across synthetic worker and rank layouts.

## Verification

```sh
make test-turbulence-core-unit
make test-reproducibility-unit
```

The first command executes registered `TURB02`–`TURB23` and `TURBOWN01` plus
CLI/source gates. `TURB21` compares a nonuniform periodic profile with exact
translated cell averages at two resolutions; `TURB22` compares the complete
time history with the analytical integral of a prescribed linear-in-time
wave-energy source. This tests source application and ledger accounting, not
the QLT calculation of that source from particles; the latter remains in the
native `growth_rate_validation_test.cpp` integration test. The
`TURB23` case then closes the weighted kinetic-energy change from controlled
outward and inward wave-frame scattering against opposite signed increments
in the resonant minus and plus wave branches. It requires both the total
particle-plus-wave energy and `particleExchangeJ` ledger to agree within
`2e-12 J`, with no positivity limiting. The
production catalog's established `TURB01` remains the independent 1-AU
magnetic-pressure closure. The second command executes PAR01–PAR05 plus
source-policy gates. Both compile the exact production
cores as C++11 with strict warnings and sanitizers. A complete AMPS checkout is
still required for native OpenMP/MPI evidence, SWMF handoff integration, and
long-run scientific conservation/observational campaigns.

The same controlled turbulence descriptors can be selected in a linked build:

```sh
../amps --test TURB21 --test TURB22 --test TURB23
../amps --test-group turbulence --test-json turbulence-results.json
```

Explicit group selection includes every registered turbulence case. Routine
`--all-tests` selection remains governed by each descriptor's runtime class.
