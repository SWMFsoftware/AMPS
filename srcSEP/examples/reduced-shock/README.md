# Reduced shock-front ambient coupling in `srcSEP`

This example couples the shared
`src/models/sep_corona_swcme/shock_front` provider to native `srcSEP`
field-line storage. It is a background-only case: the provider supplies a
time-dependent prescribed front plus undisturbed ambient plasma and IMF, and
`srcSEP` installs that ambient on every native field-line vertex. All particle
injection callbacks are disabled and the application verifies the actual
global AMPS particle count after every step.

The case is synthetic and physically plausible, not an observational event
validation. It qualifies this selected reduced profile and its application
boundary; it does **not** qualify BG3D-4, a sheath/ejecta volume, particle
acceleration, or SEP transport.

## Files and identity

- `field_line.in` is the complete schema-3 `srcSEP` input. The line follows
  the +Z HCI direction from 20 solar radii to 1 AU with 513 vertices.
- `positive_1au.event` is the strict reduced-provider event. Unknown,
  duplicate, missing, non-SI, or unsupported values fail before AMPS mesh
  allocation.
- `assets/synthetic_dipole.harmonics.csv` is the checksummed real-orthonormal
  PFSS asset. Its byte checksum and every resolved physical parameter enter
  the event fingerprint.

The currently resolved event identity is
`0c3e100d11b4c1c47d46a2ba8c5063966528c4cc59575fbd14936ad8cf08fe1b`.
Do not copy the event without its relative `assets/` directory.

## Physics and conventions

All dimensional values are SI. Coordinates and vectors use the heliocentric
inertial (HCI) Cartesian basis, the Sun is at the origin, and outward radial
and shock-normal speeds are positive away from the Sun. The example axis is
`q=(0,0,1)`.

The ambient authority is the maintained PFSS/Parker composite. Inside the
source surface it evaluates the checksummed potential field and the declared
coronal plasma branch. Outside it, the same magnetic flux is continued into a
retarded rotating Parker IMF. The selected open-wind density conserves radial
mass flux,

```text
rho(r) U_r(r) r^2 = constant,
```

and the pressure is computed from the declared electron/proton temperatures
and composition. The native line begins at 20 solar radii, so all of its
vertices sample the outer Parker branch; the front itself starts in the low
corona and is still evolved there by the shared provider.

The prescribed finite-SSE front has fixed direction and half width
`lambda=pi/4`. For apex radius `Ra`, its generating sphere has

```text
c = Ra/(1+sin(lambda)),   a = c sin(lambda),
Vn = Va (d.n + sin(lambda))/(1+sin(lambda)).
```

The apex starts at `8.04e8 m` with speed `100 km/s`. During the first 1800 s,
its speed follows the C2 quintic ramp

```text
Va(t) = V0 + (V1-V0) [10 s^3 - 15 s^4 + 6 s^5],  s=t/1800 s,
```

to `V1=1400 km/s`; the analytically integrated radius is used rather than a
time-step integral. The front then coasts to the handoff at exactly 20 solar
radii. The outer phase is matched in apex radius and speed and follows

```text
dV/dt = -Gamma (V-w)|V-w|,
V = w + (Vh-w)/(1+Gamma |Vh-w| (t-th)),
R = Rh + w dt + sign(Vh-w) log(1+Gamma |Vh-w| dt)/Gamma,
```

with `w=400 km/s` and `Gamma=5.037350032041981e-12 1/m`. The handoff is C1;
acceleration need not be continuous. At each surface face the shared provider
independently evaluates fast-mode admission and the conservative ideal-MHD
Rankine--Hugoniot jump. Those immediate one-sided states remain surface data;
this adapter never extends them through a downstream volume.

## Native storage, epochs, and failure semantics

At each 600-s background epoch the adapter first prepares the complete shared
front transactionally, then evaluates all 513 native field-line vertices into
a staging array. Native arrays are changed only after every query succeeds.
The outgoing magnetic field, velocity, number density, proton temperature and
total thermal pressure are copied into AMPS' previous-epoch datums before the
new realization is installed. Generation one copies the new state into both
slots because no earlier physical event epoch exists.

An immutable `reduced-shock`, model-owned snapshot binds the event identity,
epoch, validity interval and generation to those arrays. Publication is
forbidden during a mover read phase. Off-cadence, reversed, out-of-support, or
failed provider candidates leave the last committed provider metadata and
native arrays unchanged.

The maintained standalone build advances the authoritative application clock
through two AMPS update hooks per outer driver iteration. Therefore this deck
uses `run.time_step_s=300` so each outer iteration advances exactly one
600-s provider cadence. This is verified by contiguous generation logs; it is
not hidden by rounding an off-cadence request.

Every MPI rank owns the replicated field-line representation and evaluates
the same provider input. Before publication, an `MPI_Allreduce` requires exact
min/max agreement for the endpoint magnetic field, flow, density and pressure.
This is replicated field-line evidence, not an AMR owner/received-ghost claim.
After every native step a separate reduction requires the actual allocated
particle population to be zero on all ranks.

## Required native build

Run these commands from the AMPS root. The process listing must show that no
build or test is using `build` before removal.

```bash
cd /home/vtenishe/Mars2/AMPS
pwd
ps -C make -C gmake -C g++ -C gcc -C cc1plus -C mpiexec -C mpirun -C amps \
  -o pid=,ppid=,stat=,etime=,args=
rm -rf -- build
./Config.pl -application=test/sep_parker_spiral__field_line
make -j16 amps
```

The removal target is only `/home/vtenishe/Mars2/AMPS/build`. Do not remove
source, model build directories, site configuration, installed libraries, or
either Mars checkout.

## Run commands

Initialization-only check:

```bash
mpiexec -n 1 ./amps \
  --input srcSEP/examples/reduced-shock/field_line.in \
  --reduced-shock-event srcSEP/examples/reduced-shock/positive_1au.event \
  --initialization-only \
  --initialization-output-dir test_output/srcsep-reduced/init-rank1
```

Two-epoch native smoke test:

```bash
mpiexec -n 4 ./amps \
  --input srcSEP/examples/reduced-shock/field_line.in \
  --reduced-shock-event srcSEP/examples/reduced-shock/positive_1au.event \
  --total-iterations 2 \
  --coupling off --cascade off --reflection off --shock-injection off
```

Full low-corona-to-1-AU background campaign:

```bash
mpiexec -n 4 ./amps \
  --input srcSEP/examples/reduced-shock/field_line.in \
  --reduced-shock-event srcSEP/examples/reduced-shock/positive_1au.event \
  --total-iterations 208 \
  --coupling off --cascade off --reflection off --shock-injection off
```

Iteration 207 publishes time `124200 s`, generation 208, where the prescribed
apex is at 1 AU and the shared provider reports an accepted apex shock. The
last step advances the application clock to `124800 s`; that later time is not
the observer-arrival claim.

Expected provenance lines include:

```text
REDUCED_BACKGROUND epoch_s=1.020000e+04 generation=18 ...
REDUCED_BACKGROUND epoch_s=1.242000e+05 generation=208 apex_radius_m=1.495979e+11 ... apex_shock_accepted=1 ...
REDUCED_PARTICLES count=0 epoch_s=1.248000e+05
```

## Focused tests

The dependency-light adapter test exercises the production adapter and shared
archives without claiming native AMPS storage:

```bash
make -C srcSEP test-reduced-shock-background-unit
```

It checks strict event/asset resolution, identity, epoch/generation validity,
a finite physical ambient sample, rejection and rollback of an off-cadence
candidate, and outward advancement at the next epoch. Parser selection is
also covered by:

```bash
srcSEP/test/run_step1_tests.sh
```

The native commands above are required in addition to these portable tests.

## Implemented versus validated

| Capability | Implementation | Current evidence |
|---|---|---|
| Shared event/asset resolution and ambient queries | Implemented in the thin production adapter; physical equations remain in `src/models` | RSHSEP01--10: 10 PASS, 0 FAIL, 0 SKIP, 0 ERROR |
| Native current/previous field-line plasma and IMF | Implemented for all vertices with transactional staging and exact mover-facing readback | Final 1-/4-rank 1-AU logs below; all 208 generations pass readback |
| Coherent epochs and MPI state | Explicit 600-s provider generations and complete endpoint min/max reduction | One-/four-rank reduced traces are byte-identical, SHA-256 `0a7b52be95e25c1ea7cdafe2b6c55e9e7c64f18fb297fc40712fc4e2c7ca3166` |
| Zero particles | All native injection callbacks are conditionally disabled only in reduced mode; actual buffer count is reduced after every step | 208/208 checks pass at one rank and 208/208 at four ranks |
| Low-coronal launch, handoff, and 1-AU trajectory | Implemented by the shared provider; `srcSEP` consumes its epoch and ambient | 10,200-s handoff and 124,200-s accepted 1-AU apex appear in both native logs |
| Native `srcSEP` checkpoint/resume | Not implemented or validated in this extension | Unqualified; the shared or `srcSEP3D` restart evidence is not substituted |
| Sheath/ejecta volume and particles | Intentionally not implemented by the reduced profile | Deferred; BG3D-4 remains unqualified |

Final native evidence is retained at:

```text
test_output/srcsep-reduced/one-au-rank1-readback/mpi/1/rank.0/stdout
test_output/srcsep-reduced/one-au-rank4-readback/mpi/1/rank.0/stdout
test_output/srcsep-reduced/one-au-rank4-readback/mpi/1/rank.{0,1,2,3}/stderr
```

All four final stderr files are empty. At 1 AU the actual native endpoint
readback is `B=(0,0,1.257375e-9) T`, `U=(0,0,6.012234e5) m/s`,
`n_e=7e6 m^-3`, `T_p=1.4e6 K`, and `p=2.706072e-10 Pa`.

## Inputs and changes

- Change trajectory, front geometry, ambient normalization, numerical surface
  resolution, endpoint, or magnetic asset in `positive_1au.event`. Update the
  asset SHA-256 when its bytes change. Every such change creates a new event
  identity and needs fresh shared-provider qualification.
- Change line/domain/mesh resolution, observers or output names in
  `field_line.in`. Keep the line inside event radial coverage and keep the
  application epoch increment equal to `run.background_dt_s`.
- The `[swcme]` section is retained because schema 3 and the baseline Parker
  geometry require it. With `--reduced-shock-event`, it is not the selected
  vertex background authority and it does not supply a downstream CME volume.
- `macroparticles_per_step` remains positive because that is the unchanged
  baseline schema contract. Reduced mode disables the callbacks before native
  stepping and checks the real population; changing this input does not enable
  particles in this example.

## Outputs and limitations

Initialization-only mode writes the three configured Tecplot products under
the requested output directory. The field-line product is geometry, and the
AMPS data product is the maintained Cartesian mesh output. Runtime coupling
evidence is emitted as the `REDUCED_BACKGROUND` and `REDUCED_PARTICLES` lines;
this change does not introduce a new field-line plasma serializer.

Known limitations are deliberate:

- the volume background is undisturbed ambient everywhere; there is no sheath,
  ejecta, wake, contact, draping, or shocked-material memory;
- immediate downstream values and Mach/compression diagnostics exist only on
  the shared front epoch and are not copied into `srcSEP` vertices;
- the front is prescribed rather than dynamically driven by the plasma;
- the synthetic dipole and isothermal ambient are controlled references, not
  a forecast or observation-derived reconstruction;
- the field line starts at 20 solar radii, so `srcSEP` does not spatially
  sample the provider's low-coronal segment even though the same event evolves
  continuously from its low-coronal launch;
- particle source, mover, scattering, wave growth and feedback paths remain
  baseline code and are intentionally inactive in reduced mode;
- native checkpoint/resume of this new `srcSEP` provider has not yet been
  qualified. Shared-provider serialization or the qualified `srcSEP3D`
  restart cannot substitute for an actual `srcSEP` restart test.
