# Reduced SSE front and Parker material line in `srcSEP`

This example uses the existing core provider in
`src/models/sep_corona_swcme/shock_front`. All plasma, field, trajectory and
local shock-jump equations remain there. Stages 0–4 integrate its state,
material geometry, finite-front connectivity and the three existing transport
solvers. Stages 5–6 add a separate upstream particle consumer; their executed
qualification gates are recorded in the stage status document. The reduced
provider still has no downstream transport volume.

## Inputs and physical authority

`field_line.in` is schema 6. Its configured physical timestep stays **300 s**;
its initial 513 points span the maintained +Z line from 20 solar radii to 1 AU.
`positive_1au.event` resolves the HCI/SI core profile and its checksummed
`assets/synthetic_dipole.harmonics.csv`. Keep the event and relative assets
together. The physical event identity remains
`0c3e100d11b4c1c47d46a2ba8c5063966528c4cc59575fbd14936ad8cf08fe1b`.

The native line is traced through the core B rather than given a separate
application Parker equation. Existing vertices follow the actual core U,
including radial variation and azimuthal motion. Core density conserves
`rho U_r r^2`; its pressure and electron/proton/alpha populations use the
resolved composition and thermal definitions. Native particle density is
explicitly **proton number density**, not electron density. No assumption of
`ne=np` is made for a helium-bearing event.

The existing finite SSE front has axis `(0,0,1)` and half width pi/4. With
apex radius Ra, its generating sphere has centre `c=Ra/(1+sin(lambda))` along
the axis and radius `a=c sin(lambda)`. The core supplies normal speed
`Vn=Va (axis dot n+sin(lambda))/(1+sin(lambda))`. A generating-sphere rear root
is not a shock. The consumer validates the leading/support branch and obtains
all local Mach, compression and jump states through the core evaluator.

The existing analytical history starts at 8.04e8 m and 100 km/s, follows a
1800-s C2 quintic speed ramp to 1400 km/s, then matches the configured outer
drag history at 20 solar radii. This integration does not alter that history.
Immediate upstream/downstream RH limits remain front data; no sheath, ejecta,
wake, shocked-material history or downstream volume is supplied.

## Time, geometry and sampling

One AMPS step advances the physical clock once by 300 s. The core's separate
prepared-front cache cadence is 600 s. Immutable explicit-time queries consume
the actual full physical interval, including a cadence crossing. They do not
round the particle clock, shorten a step or add a second clock increment.
Old documentation describing two clock increments per iteration is superseded.

After particle reads join, the native owner prepares every moved vertex,
its core state, target connectivity and publication metadata. A new vertex's
previous state is evaluated on the inverse material flow, not copied from a
neighbour. Core-cache candidates and metadata are preallocated before MPI
readiness; failures do not publish a partial native state. Replicated physical
states must agree exactly across ranks before commit.

Every endpoint moves. Inner extension follows B by a full first-segment arc
increment and repeats until the minimum radius is covered. Tracing along -U
would generally leave a Parker magnetic field line. Inserted segments start
empty with zero sampling exposure. An entirely exterior tail chord is removed
with its particles; represented losses and retired local sampling history are
recorded. Per-segment exposure preserves the unchanged global output schedule
without treating a new segment as if it existed for the whole sample window.

## Build and bounded run

Use an isolated configured CPU field-line checkout; preserve any active build.
The tested configuration is:

```sh
./Config.pl -spice-path=nospice -spice-kernels=nospice \
  -application=test/sep_parker_spiral__field_line-fte
make -j4 share_lib
make -j4 amps
```

The retained input contains a self-consistent turbulence selection. That
unfinished polarized-production closure is still rejected. Stages 0–4 use
explicit **prescribed**, read-only turbulence/background and disable local
wave evolution/feedback. From the configured checkout, with a fresh working
and output directory, run:

```sh
mpiexec -n 4 /path/to/configured/amps \
  --input /path/to/AMPS/srcSEP/examples/reduced-shock/field_line.in \
  --reduced-shock-event /path/to/AMPS/srcSEP/examples/reduced-shock/positive_1au.event \
  --total-iterations 4 --no-test-manager \
  --turbulence-source prescribed --turbulence-coupling-policy disabled \
  --coupling off --cascade off --reflection off \
  --turbulence-advection off --shock-injection off
```

Expected geometry epochs are 300, 600, 900 and 1200 s; `REDUCED_PARTICLES`
must report the actual allocated global population as zero at each endpoint.
`REDUCED_GEOMETRY` reports insertion/removal and material generation. Retired
sampling archives are per rank: exposure is elapsed sampling opportunity,
not a quantity to add across replicated geometry ranks.

A 1-AU background campaign would now need 416 complete 300-s steps to reach
124800 s; the prescribed apex reaches 1 AU at 124200 s. Historical 208-step
receipts used the old clock behavior and are not qualification of this moving
owner. The stage-0–4 qualification is bounded; it does not claim a new complete
1-AU, coupling, restart or observational campaign.

## Validation commands and scope

```sh
make -C srcSEP test-reduced-shock-background-unit
make -C srcSEP test-standalone-integration-unit
python3 srcSEP/test/run_standalone_integration_tests.py \
  --sanitizers --output-dir /tmp/sep-integration-FRESH
python3 srcSEP/test/run_native_standalone_integration.py \
  --build-root /path/to/matching/configured/AMPS \
  --output-dir /tmp/sep-native-FRESH --ranks 1 4
```

The native fixture links the ordinary configured AMPS archive and uses actual
MPI, PIC particles, field-line pools, sampling and collective readiness. Its
explicit manufactured initial particles test transport and loss accounting;
they do not select a production injection law. Pure tests additionally grade
composition, ideal induction/continuity, curved-line refinement, finite-front
root support, tube-area drift and common diffusion limits against independent
references. See [the stage report](../../validation/STANDALONE_INTEGRATION_STATUS.md)
and [test contract](../../validation/STANDALONE_INTEGRATION_TEST_CONTRACT.md).

## Remaining limits

The default reduced mode retains the zero-particle invariant. The core event's
`particle_mode=disabled` means the core does not run particles; an explicitly
enabled application-owned consumer has its own fingerprinted contract.
Downstream continuation, PFSS material evolution, moving wave remap,
nonzero feedback, moving-topology checkpoint/restart and other core backends
are not admitted by these tests. Existing SWCME schema fields remain required
input compatibility; `--reduced-shock-event` selects the actual core authority.
No reference solution or physical gate is regenerated to accommodate results.

## Explicit upstream particle consumer

Compile the source-enabled example with the application-specific proton
profile. The legacy `H_PLUS` species catalog uses 1 amu, which is not the
proton mass. This wrapper first retains that catalog's inactive-species macros
and then overrides only the modeled species table with the physical proton
mass. It does not change the shared species catalog or other applications.

```sh
./Config.pl -spice-path=nospice -spice-kernels=nospice \
  -application=test/sep_reduced_front__field_line
make -j16 share_lib
make -j16 amps
```

Add the following options to the prescribed, read-only command above:

```sh
--reduced-shock-particles on --shock-release-distance-m 1e8 \
--shock-return-policy delete --shock-source-normalization net-upstream-escape \
--total-iterations 48
```

Here **1e8 m is an explicitly manufactured verification distance**, not a
calibrated production default. The profile's 300-s physical timestep and
400-step output sampling window remain unchanged. Forty-eight steps exercise
the beginning of the connected source; a complete sampling-output window
requires its original 400 physical steps. Holding L_ref fixed during numerical
refinement is mandatory. Record lambda/L_ref and inflow/particle-speed scales
when judging the physical applicability of a Parker diffusion approximation;
numerical transport tests do not establish that approximation for every input.

The `source: shock` energy bounds are total proton kinetic energies in joules.
For this consumer `injection_efficiency` prescribes the net first-passage
number fraction eta, while `samples_per_step` sets stratified temporal proposal
count N. The core local compression gives q=3r/(r-1), and the reused source law
samples **dN/dp proportional to p^(2-q)**. No SWCME ambient density, compression,
or total-front area is borrowed through its compatibility source helper.
Other required SWCME schema fields remain compatibility inputs; the reduced
core and compiled proton species are the physical authorities in this mode.

For conserved tube flux Phi, the connected shock patch is Phi/|B dot n|.
The prescribed net number rate is eta*n_p*(V_n-U_1 dot n)*Phi/|B dot n|.
It is a rate for this tube in physical protons per second. Source energy is
the explicitly injected already accelerated population; the application does
not resolve accelerator backreaction or create an energetically modified RH
jump. The transport ledger treats that source kinetic energy as an external
input, separately from later deterministic and wave-frame work.

Birth j uses an independent keyed open uniform in its temporal stratum and
weight R(t_j)*dt/N. Actual-epoch source geometry intersects the same moving
line with the upstream positive normal-distance offset. Focused pitch angles
sample positive normal flux from an upstream-plasma-frame isotropic parent;
shock-frame isotropy is not substituted. Parker integrates that angular
source, with an isotropic diagnostic projection. Configured `1e8` remains a
normal distance, not an along-line displacement or a grid-dependent epsilon.

Sub-fast or disconnected proposals create no particles. No-focused-escape
exclusions are recorded separately, without redistributing their weights.
Returning particles are absorbed at the first resolved encounter and retain
distinct birth and event-energy accounting. `pass` and `gross-emission` are
recognized but terminate with explicit unsupported messages. No missing
downstream or gross-to-net escape prescription is invented.

`REDUCED_SOURCE` records candidate, excluded and committed release measures.
`REDUCED_SOURCE_LEDGER` records retained, returned, domain-exited and retired
particles, their birth/event energies, signed work and conservation residuals.
Retained birth energy is printed separately from current kinetic energy so an
external checker can audit both balances without assuming that cooling or
prescribed scattering conserves particle kinetic energy by itself.
PIC samples the departure population before its source callback, so future
births cannot enter an earlier sample. This ordering uses the frozen birth
window; no physical clock or sampling denominator is adjusted to fit it.

The active solution is exported in `amps.FieldLines.out=*.dat`, its material
`.sampling` sidecar and the reduced source/energy ledgers. The separately named
`swcme-comparison-shock-vs-time.dat` is the compatibility SWCME curve, rather
than the selected reduced front. `background.fl=*.dat` tabulates the explicitly
named alternative mean-free-path formulas for comparison; those formulas do
not silently replace the configured production coefficient provider.

The first consumer supports one transverse leading barrier. Tangency,
simultaneous multiple absorbing fronts, occupied connectivity/criticality
changes and activation of a merger without the birth-energy invariant are
explicit unsupported cases. They are not silently passed through. Stage 7
retains distributed migration, cohort-preserving resampling and moving restart
qualification. Use `--reduced-source-verification-energy-j` only for an
explicitly labelled manufactured monoenergetic test; zero selects core DSA.
