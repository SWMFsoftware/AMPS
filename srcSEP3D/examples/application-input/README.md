# Shared-input reduced-front example

This directory is a self-contained example of the one-dash srcSEP3D runtime
input. `amps.in` is the global container, `sep3d.in` owns application controls,
and `reduced-shock-surface.in` contains the complete strict reduced model. Its
PFSS asset is the checksummed synthetic dipole already maintained under
`../shock-front/assets/`; there is no external data prerequisite.

The case uses SI units and HCI coordinates. It launches a 45-degree half-width
finite SSE front at `8.04e8 m`, accelerates smoothly from `100 km/s` to
`1400 km/s` during 1800 s, hands continuously to the quadratic-drag branch at
`20 R_sun`, and declares an accepted-shock observer at 1 AU. The ambient is the
maintained PFSS/Parker reference with a checksummed dipole, 7 cm^-3 electron
density at 1 AU, 1.4 MK species temperatures and no alpha abundance. This is a
synthetic plausibility case, not observational validation.

`amps-injection-smoke.in` is a separate native allocation fixture. It starts
with an already-formed, constant-speed 1400 km/s front at 50 R_sun, hands over
at 51 R_sun, and runs exactly two steps. The late start keeps its complete SSE
cap outside the application's 20-R_sun particle-source shell. That artificial
initial condition is useful for proving accepted-face sampling and MPI
ownership without waiting through formation, but it is not a coronal launch
model and must not replace `amps.in` in physical campaigns.

Build from the AMPS root (the root `build` deletion is mandatory before every
native rebuild):

```bash
pwd -P
pgrep -af '(^|/)(make|g\+\+|mpiexec|amps|stage1|run_tests\.py)( |$)' || true
rm -rf -- build
./Config.pl -application=sep3d
./ampsConfig.pl -input sep3d.input -no-compile
make -C srcSEP3D prepare-production
make -j16 amps
```

Run a complete native initialization on one rank:

```bash
mpiexec -n 1 ./amps \
  -input srcSEP3D/examples/application-input/amps.in \
  --initialization-only \
  --initialization-output-dir test_output/application-input-model/native-rank1
```

Use `mpiexec -n 4` and a distinct output directory for the MPI agreement run.
With no CLI arguments the program instead reads `./amps.in`.

The checked-in example has `maximum_time_steps=4`, making a direct invocation
a bounded smoke run. Its launch begins sub-fast, so these first four steps
exercise the valid zero-rate Poisson state and intentionally inject no
particles. Increase the horizon to cover physical shock formation; do not add
a Mach/compression floor merely to manufacture startup events.

Exercise a positive native source immediately with:

```bash
mpiexec -n 1 ./amps \
  -input srcSEP3D/examples/application-input/amps-injection-smoke.in \
  --output-dir test_output/srcsep3d-shock-injection/positive-rank1

mpiexec -n 4 ./amps \
  -input srcSEP3D/examples/application-input/amps-injection-smoke.in \
  --output-dir test_output/srcsep3d-shock-injection/positive-rank4
```

The rank-zero source lines must agree in event fingerprint, step, generation,
accepted-face count, total physical rate, Poisson candidate count, allocated
all-rank count, and disconnected count. A zero disconnected count plus
`injected_all_ranks=poisson_candidates>0` proves that every sampled point was
owned and allocated exactly once; it does not validate an acceleration
efficiency or the artificial prehistory of this smoke fixture.

The source preflight remains reproducible across one and four ranks. Both
layouts produced the same immutable normalization receipt:

```text
h_min                         = 1.34953704560687232e9 m
dt                            = 2.02430556841030818 s
normalization epoch           = 10199.9999999989086 s
accepted curved shock area    = 2.67211070324907966e20 m2
electron incident rate        = 3.25066729609523584e35 s-1
electron macroparticle weight = 6.58034390853486657e32
configuration fingerprint     = c5ffd8a154b2277f
```

That generated `sep3d` build contains one compiled species, `ELECTRON`.
Proton and alpha weights are computed independently when those populations
are present in both the generated AMPS species table and the ambient
composition; an absent/zero population is a configuration error rather than a
borrowed electron or proton weight. Native logs and initialization products
are under
`test_output/application-input-model/native-rank{1,4}-final-v2/`.

Run the dependency-light parser/numerics and shared-physics regressions with:

```bash
make -C srcSEP3D test-stage1 -j16
make -C src/models/sep_corona_swcme -j16 shock-front-test
python3 srcSEP3D/test/run_tests.py --all --amps-source . --no-build \
  --output-dir test_output/srcsep3d-shock-injection/aggregate-with-source-final
```

The post-source results are respectively `PASS=151 FAIL=0 SKIP=0 ERROR=0`
and `PASS=57 FAIL=0 SKIP=0 ERROR=0`; the complete application aggregate is
`PASS=173 FAIL=0 SKIP=29 ERROR=0`. Its 29 SKIPs are declared external/native/
observational evidence gates, not failures of this source. Current native
one-/four-rank evidence is recorded in `CODEX_REDUCED_SHOCK_PLAN.md`; the
earlier initialization-only runs do not qualify native injection.

Edit `sep3d.in` to change the finite iteration horizon, target mean
model-particle count at the
normalization radius, kinetic-energy interval, phase-space power law, fatal
Poisson guard, maximum represented speed, CFL margin or required
source-normalization radius. Edit the reduced subsection
to change launch, handoff, ambient or magnetic parameters. Any magnetic asset
change requires updating `assets.harmonics_sha256`; an incorrect checksum is a
fatal data-integrity error. The normalization radius must lie on the declared
monotonically outward trajectory and contain a positive accepted-shock area.

Startup reports `dt=f h_min/v_max` and, for each compiled species,
`W_s=Ndot_s dt/N_model`. The gross `Ndot_s` is the accepted upstream incident
flux and no hidden acceleration efficiency is applied. At each step the AMPS
user-defined injection callback recomputes accepted triangular-face rates and
draws exponential waiting times at `lambda=Ndot_live/W_s`. Therefore
`particles_per_iteration` is the expected count only at the normalization
radius, not a forced count at every epoch. Particles are sampled uniformly on
the selected planar triangle, launched anti-sunward, and allocated only by the
owning MPI rank. Rank zero prints the total live rate, Poisson candidate count,
all-rank allocated count and disconnected-domain count per species and step.

`RequestedParticleBufferLength=1500000` in the generated application input is
capacity allocated independently on every MPI rank; it is not a source-rate
limit and does not replace the fatal per-step Poisson guard. In the two-step
positive fixture the observed peak memory was about 1.00 GiB for one rank and
2.61 GiB total for four ranks (about 0.65 GiB per rank). Reduce this capacity
only after bounding the largest supported burst, because AMPS treats particle
buffer exhaustion as a fatal allocation error.

Production preparation installs the srcSEP3D mover as the *last effective*
`_PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_` definition in the
generated AMPS configuration. This ordering matters: AMPS appends a generic
guiding-centre default after the application template, and an earlier macro
definition would leave injected particles invisible to the srcSEP3D source
ledger. `make -C srcSEP3D audit-production-symbols` and `BLDL3D05` check the
effective final definition rather than accepting a stale earlier definition.

This is a synthetic seed source based on gross swept-up ambient population; it
is not an observationally calibrated acceleration efficiency. The provider
remains ambient-only away from the surface: no output value should be
interpreted as a sheath or ejecta plasma state.
