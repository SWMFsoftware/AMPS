# Configured AMPS core-test application

This is the native test boundary for core behavior that cannot be compiled as
a dependency-light unit test. It deliberately uses the small harness in
`src/test/native/test_harness.h`, so no machine-specific GoogleTest include or
library path is required.

The application is not a physical Enceladus, exosphere, or plasma model. Its
input selects the seven neutral species needed by the H2O photolysis and
internal-energy fixtures. The target-free, SPICE-free `exosphere.pl` block
installs AMPS' default internal-boundary user-data type definitions required by
mesh headers; all exosphere sources and the exosphere module remain excluded.
Most cases are process-local. The physics fixture lazily builds a small 3-D
AMR mesh with eight allocated owner blocks, registers one analytic sphere,
allocates real particles in mesh cells, and invokes production movers, NTC
collisions, Larsen-Borgnakke relaxation, and photolytic decay scheduling.

Do not build it manually while another AMPS build or run is active. The global
runner performs the repository preflight, removes only the generated `build`
directory, selects this application through `Config.pl`, builds with `-j16`,
and preserves structured evidence:

```sh
python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001
python3 src/test/run_tests.py --case PIC-UNITS-GLOBAL-INIT-001 \
  --case PIC-UNITS-CONVERSION-TABLE-001
python3 src/test/run_tests.py --test CORE-PIC-WEIGHT-TIME-001
python3 src/test/run_tests.py --campaign profile
python3 src/test/run_tests.py --test CORE-PIC-SPLIT-MERGE-001
```

The current cases verify process-local particle-buffer allocation/deletion,
state access and cloning, and intrusive-list operations. Rank zero executes
these local cases while all ranks participate in MPI initialization and final
status reduction. Future ownership/halo tests must be separate collective
cases that inspect owner and ghost state on every rank; rank-zero execution is
not sufficient for those tests.

`ParticleBufferScope` saves every process-global particle-buffer pointer and
counter, replaces them with a 64-record fixture, releases exactly the storage
allocated by that fixture, and restores the original state. On the Linux
non-CUDA path AMPS uses `aligned_alloc`; the fixture correspondingly uses
`free`, not `delete[]`. This ownership distinction is part of the test design,
because a mismatched deallocator can make a passing assertion corrupt later
cases.

The final rank status is an MPI maximum over the harness vocabulary
`0=pass`, `1=assertion failure`, `2=usage/no-selection error`, so every rank
leaves MPI with the same dominant outcome.

## T2 particle records, weight/time profiles, and split/merge

T2 is qualified for six one-rank native profiles and 85 profile cases. The
baseline `CORE-PIC-NATIVE-001` owner contains 43 cases, including seven packed
fixed-layout cases and five species-global weight/time cases. The five other
owners are `CORE-PIC-SPLIT-MERGE-001` (13 packed individual-weight cases),
`CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001` (7),
`CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001` (13),
`CORE-PIC-WEIGHT-TIME-LOCAL-001` (5), and
`CORE-PIC-WEIGHT-TIME-SINGLE-GLOBAL-001` (4).

`particle_layout_t2.cpp` is compiled in four record layouts. The packed
profiles use the production optional fields; the aligned profiles additionally
request anonymous 1-, 2-, 4-, 8-, and 16-byte fields with matching alignments.
Cases independently calculate each aligned offset and record stride, then
check clone, delete/reuse with payload clearing, intrusive-list integrity,
capacity, invalid host access, and real exhaustion. Exhaustion runs in its own
MPI process and passes only when the process exits nonzero and reports
`The particle buffer is full`. The production request path validates power-of-
two alignments and arithmetic overflow, and `SetParticleDeleted` addresses the
requested record instead of the buffer base.

The weight/time profiles cover species-global, species-local block-associated,
and single-global storage. They check round trips, minimum reduction where
applicable, scaled cross-species distribution copies, local physical-
distribution rescaling, shared-slot aliasing, and typed host rejection of
invalid species indices and zero, negative, NaN, or infinite values. Copy
operations reject the `-1` uninitialized source sentinel and finite overflow
before modifying the target. Individual particle-weight accessors apply the
same finite-positive host contract. Public block accessors remain host/device
functions, but CUDA/device invalid-state diagnostics are not T2-qualified.

`split_merge_t2.cpp` ports the fixed- and variable-weight coverage to both
packed and aligned individual-weight profiles. Independent before/after
moments cover model count, physical weight, species, weighted position,
momentum, translational energy, rotational energy, and every vibrational mode,
with empty/one/two-particle, target-unreachable, boundary, and unrelated-
species bins. Generic code cannot infer the semantics of anonymous application
bytes, so split uses an exact source clone and merge retains an exact survivor;
the tests accept either valid source record rather than inventing an average.

The populated-cell fixture may observe AMPS's `-1` pre-initialization
sentinel.  Its transactional teardown restores that raw saved global or local
storage directly instead of passing lifecycle state through public physical-
value setters, whose minimum-reducer behavior cannot perform an exact restore.
This distinction was found by the first T2 regression campaign and is covered
by the subsequent baseline and six-profile order-isolation receipts.

The qualification receipts are `/tmp/amps-t2-native-full-normal` (43/43),
`/tmp/amps-t2-profile-reverse` (85/85), and
`/tmp/amps-t2-profile-shuffle-1729` (85/85). Focused final receipts are
`/tmp/amps-t2-layout-packed-fixed-recheck` (7/7),
`/tmp/amps-t2-split-merge-packed-recheck` (13/13),
`/tmp/amps-t2-layout-aligned-fixed-initial` (7/7),
`/tmp/amps-t2-layout-aligned-individual-recheck` (13/13),
`/tmp/amps-t2-weight-local-recheck` (5/5), and
`/tmp/amps-t2-weight-single-global-initial` (4/4). Stage T3 mesh/AMR work
remains planned; T2 does not qualify restart or distributed ownership.

Immediately before each selected C++ test function, the shared harness prints
and flushes `[CASE] RUN  <ID> <description>`. Both the native adapter and global
runner stream that record, so a long statistical case identifies itself before
execution. Final `[ID] PASS|FAIL` records remain the structured evidence parsed
by the adapter.

## T4 movers, boundaries, injection, and sampling

The baseline executable owns 19 T4 cases and a separate periodic compile
profile owns one. Constant SI fields are supplied through configured hooks;
the production Boris pushers are compared with analytic Lorentz solutions at
multiple step sizes. The relativistic case uses constant gamma and
`omega=qB/(gamma*m)`. The guiding-center test calls a public field-local RHS
and independently evaluates parallel electric/mirror acceleration plus
`E x B`, grad-B, and curvature drifts. Neutral species, zero magnetic field,
negative magnetic moment, and non-finite input return typed statuses and zero
outputs instead of plausible-looking drift values.

Boundary fixtures solve their event geometry independently. Static sphere
cases cover specular reflection, absorption, tangent/no-hit, start-inside, and
outer loss. Moving plane and sphere cases use an inertial frame, solve the
linear or quadratic relative-motion event exactly, then reflect relative to
the local surface velocity. A separate clean periodic build verifies real/
ghost pairing, wrapping, intrusive-list transfer, and associated-node-data
synchronization.

Injection checks use the half-space Maxwell flux, exponential waiting-time
mean/variance, uniform Cartesian cell and barycentric triangle positions, and
the spherical-volume CDF `P(r<R)=(r/R)^3`. Sampling reconstructs raw weighted
number, velocity, flux, full pressure tensor, kinetic energy, speed, and
temperature from a literal population, then repeats with a legal particle-
count/statistical-weight rescaling. Collect/swap/flush proves that no stale
buffer state survives an interval boundary. These are one-rank contracts;
coarse/fine motion, production non-spherical/diffuse surfaces, source
ownership/rate normalization, histograms, and MPI migration remain open.

## T5 collisions and internal degrees of freedom

Twenty baseline selector cases drive NTC, simplified NTC, MF, simplified MF,
and Yinsi MF through one independent populated-cell rate oracle. Each selector
has unlike/same directed-rate, momentum/total-energy, zero-relative-speed,
empty/singleton, and aggressive collision-limiting evidence. NTC, MF, and
Yinsi also exercise unequal global weight/time-step update probabilities.
Only NTC and MF advertise and test individual record weights; simplified
selectors explicitly require fixed record weight and Yinsi reports the mode
unimplemented.

Clean VHS and VSS profiles reconstruct `sigma(g)` from literal deck values and
integrate `<sigma*g>` against the Maxwell relative-speed density without using
production expected values. VSS additionally checks
`cos(chi)=2*U^(1/alpha)-1`, its first two analytical moments, and elastic
momentum/energy. The binary alpha closure is the arithmetic mean, matching the
declared diameter mixture; it is documented as a model choice rather than a
transport fit.

LB cases cover storage, equilibrium, microcanonical RT partitioning in both
hot/cold directions, every real/virtual record-update flag, and monatomic,
linear, nonlinear, and multi-vibrational-mode limits. Three generated profiles
measure event probabilities `1/Zrot` at `Zrot=2,5,10`. A monatomic partner uses
the rotating species' value; unequal values for two rotating species are
typed ambiguous rather than silently mixed. RT is supported, VT remains
compiled but unverified and terminating, and VV is not implemented.

Final retained evidence is 51/51 in normal, reverse, and seed-1729 shuffled
baseline order at `/tmp/amps-t45-baseline-final-{normal,reverse,shuffle-1729}`
plus 7/7 selected variant cases across six clean builds at
`/tmp/amps-t45-profile-variants-final`. The complete final-tree T2-T5 profile
campaign at `/tmp/amps-t2-t5-profile-qualified` passes 194/194 cases in 12/12
suites.

## T6 photochemistry and compiled reaction models

The baseline photochemistry owner executes six cases against populated local
blocks. The fixture deliberately rebuilds the AMR local-block table after
allocating its eight owner blocks; the production chemistry driver iterates
that table, so omitting this step would turn a decay test into a zero-block
no-op. The cases compare observed H2O parent survival with the exact
`exp(-dt/tau)` law, check disabled and invalid-lifetime paths, verify both
intrusive-list links after deletions, and repeat initialization to prove that
the lifetime table is not accumulated twice. Stable `expm1`/`log1p` forms are
used by production for small probabilities and Poisson event times.

A separate clean individual-weight profile enables one daughter-producing H2O
channel. Its oracle is independent of the production channel sampler and
closes daughter identity/count, statistical weight, charge, parent/daughter
mass, centre-of-mass momentum, and excess kinetic energy. This qualifies that
configured channel only; it does not claim MPI/thread-safe insertion, a
general runtime network, or application-specific multi-parent callbacks.

The Burger reaction profile exposes three compiled parents and tests total
electron-impact rates at literal energies, exact threshold behavior,
conditioned eligible-channel selection, H2O product statistics, and typed
no-reaction handling. Lindsay--Stebbings and Maher--Tinsley tests independently
evaluate their H/H+ analytic lifetime formulas at valid points and require
the established large-lifetime sentinel outside the physical domain. These
charge-exchange files do not expose event/product/particle-weight APIs, so T6
does not infer those capabilities.

The dependency-light `CORE-CHEMICAL-SOLVER-001` suite is built outside this
configured application. It checks irreversible and reversible closed forms,
a stiff decay, a coupled chain, positivity/conservation, and typed invalid
networks using independent analytic solutions.

Run the complete bounded T6 gate with:

```sh
python3 src/test/run_tests.py --test CORE-CHEMICAL-SOLVER-001 \
  --test CORE-PHOTOCHEM-EXECUTION-001 \
  --test CORE-PHOTOCHEM-DAUGHTER-001 \
  --test CORE-CHARGE-EXCHANGE-001 \
  --test CORE-ELECTRON-IMPACT-001
```

The retained normal-build and order-isolation evidence is listed in
`src/test/STATUS.md`. The manifest aliases charge exchange and electron impact
to one exclusive reaction-profile owner, so selecting both does not rebuild
the same executable twice.

## Physics suites and independent oracles

The original `CORE-PIC-PHYSICS-001` foundation contains 19 native cases in six
groups; its compatibility alias now also selects the T4 physics extensions.
T5 has dedicated collision/IDF aliases and compile-profile owners. Every
statistical case sets a fixed Mersenne-Twister seed locally so test ordering
does not change its sample stream.

- Photochemistry kernels check the seven active-Sun Huebner H2O rates against a
  hand-entered sum, the analytical inverse-square rate at 2 AU, and the OH
  yield from channels 0 and 6. The generic decay scheduler is compared with
  `1-exp(-dt/tau)` using 100,000 Bernoulli trials and a six-sigma binomial
  bound. Two- and three-body products must close centre-of-mass momentum and
  excess kinetic energy to normalized `2e-12`; 120,000 two-body directions
  must have zero component means and component second moments `1/3` within
  `0.008`. A 240,000-draw production channel-selector test reproduces all
  seven tabulated branching probabilities within six binomial standard
  deviations. A separate product sample checks the neutral-only build's
  compact available-daughter list, including deliberately unavailable ion and
  electron products. A NaN-sentinel sparse-product case requires velocity
  vectors to be compacted into the same returned indices as daughter IDs.
- Hard-sphere collision kernels check total momentum, translational
  kinetic energy, and relative speed for unequal masses in both production
  rotation branches. A separate 120,000-event moment test requires a uniform
  direction on the sphere.
- Distribution sampling checks a drifting Maxwellian's three component means,
  the independent mean thermal speed `sqrt(8*k*T/(pi*m))`, and recovered
  kinetic temperature over 160,000 draws.
- Internal degrees of freedom check the configured particle byte offsets and
  three H2O vibrational modes, rotational equipartition
  `<E_rot>=f*k*T/2`, the harmonic-oscillator mean
  `<E_vib>=k*Theta/(exp(Theta/T)-1)`, and microcanonical beta-partition mean
  `f_internal/(f_internal+f_translational)`. Every partition draw must also
  conserve its fixed total energy. RT cases independently check the two-stage
  beta law, exact energy conservation, molecular/monatomic branch, and the
  current closure's analytic transition probability `1/Z_rot=1/5`. VT and VV
  are not executed because the current production VT path explicitly declares
  itself unverified; enabling it would turn a known unsupported path into
  false coverage.
- Movers use the exact neutral, zero-field, target-free solution
  `x(t)=x0+v*t`, `v(t)=v0`. Both production movers must reproduce it to
  `3e-15`, locate the destination AMR leaf, and place the particle at the head
  of that cell's temporary intrusive list. A third case launches a particle
  toward a radius-0.2 m sphere, checks the analytic hit time, dispatches one
  specular-reflection callback, and verifies normal-velocity reversal plus the
  remaining flight interval and final cell insertion.
- Unit normalization uses the production global `PIC::Units::Factors` and
  public wrappers in `src/pic/units/units.cpp`. The initialization case clears
  that global transactionally, initializes the AMU/e system, and checks its
  mass and charge scales against independent SI references. It then repeats
  initialization and requires exact idempotence. The print case writes through
  `PIC::Units::PrintConversionTable` to a temporary stream and verifies both
  the schema and a formatted value from the active global object. Both cases
  restore the prior global factors so normal, reverse, and shuffled execution
  have the same initial state. These cases qualify serialized startup only;
  concurrent normalization reconfiguration is not a supported API.
- Three populated-cell cases exercise combined production modules rather than
  isolated kernels. A 160-particle counterstreaming H2O population undergoes
  NTC collisions and RT equilibration; total momentum and translational plus
  rotational energy are conserved, rotation becomes nonzero, and its
  equilibrium fraction must be within `0.10` of one half. A 240 H2O plus 240 O
  Maxwellian cell checks the sampled cross-collision rate against
  `n_H2O*n_O*sigma*sqrt(8*k*T/(pi*mu))` within `0.12`. A 2,400-particle H2O
  cell applies the production scheduler and parent removal, reproducing
  `exp(-dt/tau)` survival within six binomial standard deviations while
  deliberately producing no daughters.

The original populated collision-rate tolerance includes finite-particle
temperature fluctuations, stochastic NTC acceptance, and sampled initial
Maxwellians; it is not a numerical-conservation tolerance. T5 extends that
foundation selector-by-selector as documented above. T4 extends the original
neutral/fixed-sphere mover evidence with charged fields, moving/absorbing/
periodic boundaries, injection, and physical sampling. The explicitly listed
T4/T5 gaps remain separate required work.

These cases are explicitly one-rank. Their manifest entry overrides a machine
profile's general MPI rank count, because a rank-zero-only local invariant is
not evidence for distributed ownership. Future MPI suites must inspect owner
and received-ghost state on every rank.

The two-body photoproduct angular test found a production defect: polar angle
was sampled uniformly, giving the non-isotropic moment `<n_z^2>=1/2`.
`PhotolyticReactions_Huebner1992ASS.cpp` now samples `cos(theta)` uniformly on
`[-1,1]`, which gives the solid-angle measure and the required `1/3` moments.

Add application test translation units to `MAINOBJ` in `makefile`, register
cases with `AMPS_TEST_CASE`, and add their stable case IDs to the individual
runner's expected registry. A missing case then fails the suite rather than
silently reducing nightly coverage.

AMPS places `MAINOBJ` members in the static archive `main.a`. Linkers do not
load an archive member merely for its static registration constructors. Each
new test translation unit must expose a link-anchor function and `main.cpp`
must call it; `LinkParticleBufferTests()` is the reference pattern. The call
has no runtime behavior, but it forces the object and its registrations into
the executable on GNU and other archive linkers.

Particle-buffer output contains exactly these stable records:

```text
[PIC-PBUFFER-INIT-001] PASS ok
[PIC-PBUFFER-STATE-001] PASS ok
[PIC-PBUFFER-LIST-001] PASS ok
```

MPI banners and diagnostics are retained in
`native-particle-buffer-output.txt`, but only records matching the stable case
grammar enter structured case evidence. A missing, duplicated, or renamed
record fails `PIC-PBUFFER-REGISTRY`; order may vary by compiler/linker.

Physics output must contain exactly the 19 `PHYS-*` records registered in
`src/test/tests/pic/native_core/run_test.py`. A missing, duplicated, or renamed
record fails `PHYS-REGISTRY`; order may vary by compiler/linker, and the
complete console transcript is retained as `native-physics-output.txt`.

The detailed mapping from every legacy `test/srcUnitTests` source to its
replacement, including the T2 individual-weight split/merge profiles and
identified legacy defects, is maintained in
`src/test/SRC_UNIT_TESTS_MIGRATION.md`.
