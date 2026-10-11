# Reduced-front injection example (unified runtime input)

This directory is a self-contained srcSEP3D example written in the generic
AMPS runtime input format (`#begin`/`#end` blocks, `!` comments, `\`
continuation, `#include` relative to the including file; normative
specification: [`../../../src/INPUT_FILE.md`](../../../src/INPUT_FILE.md)).
The files are:

| File | Contents |
|---|---|
| `amps.in` | complete production-like run: `run` block (derived time step, four-step horizon, nested `restart`) and `source: shock` block (accepted-shock incident-flux normalization), then `#include "sep3d.in"` and `#include "reduced-shock-surface.in"` |
| `amps-injection-smoke.in` | the same structure for the two-step native allocation fixture; includes `sep3d.in` and `reduced-shock-surface-injection-smoke.in` |
| `sep3d.in` | srcSEP3D blocks shared by both run files: `application: sep` (shared SEP schema 6, `intent = shock-injection`, `transport = parker`, nested `transport`/`population_control`/`storage`/`turbulence`), `mesh` (with `domain`, `field_line`, both `refinement` blocks, `active_region`, `memory`), `physics: background` (`provider = runtime-model`, `model_id = sep-corona-swcme-shock-front-v1`, nested `parker` with the seven ambient keys, required because there is no `physics: swcme`), `physics: shock` (`authority = none`), `observer: default`, `output` |
| `reduced-shock-surface.in` | the complete strict reduced shock-front plus ambient model as one `physics: shock_front` block |
| `reduced-shock-surface-injection-smoke.in` | the smoke fixture's `physics: shock_front` block |

The former `sep3d-injection-smoke.in` was deleted; the smoke fixture's run and
source differences now live in `amps-injection-smoke.in`.

Provenance of `sep3d.in`.  Before the input unification this example used the
retired one-dash `#section begin: sep3d` container, which set only the
shock/source keys and silently took every other value from compiled
`RunConfiguration3DOptions` defaults.  Those values are now written explicitly
with the same numbers; comments marked `[former default]` identify them.  They
are plausibility settings for this example, not tuned production choices.

The reduced model's PFSS asset is the checksummed synthetic dipole already
maintained under `../shock-front/assets/`; there is no external data
prerequisite.  A relative asset path in a `physics: shock_front` block
resolves from the directory of the file that holds the block.

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
./Config.pl -application=sep3d   # deck selects runtime mover dispatch
make -j
```

Run a complete native initialization on one rank:

```bash
mpiexec -n 1 ./amps \
  -input srcSEP3D/examples/application-input/amps.in \
  --initialization-only \
  --initialization-output-dir test_output/application-input-model/native-rank1
```

Use `mpiexec -n 4` and a distinct output directory for the MPI agreement run.
`--input` is an exact synonym of `-input`.  With no input argument the
program reads `./amps.in` in the working directory (as a coupled SWMF run
does), so copying this directory and running `./amps` inside it is
equivalent; the relative `#include` lines resolve from the including file.
`--dry-run` and `--restart` are accepted for this input like for every
other srcSEP3D input.

The checked-in example has `run.maximum_time_steps = 4`, making a direct invocation
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

Historical receipt (old one-dash format).  Before the input unification,
the equivalent one-dash deck produced the following identical normalization
receipt on one and four ranks:

```text
h_min                         = 1.34953704560687232e9 m
dt                            = 2.02430556841030818 s
normalization epoch           = 10199.9999999989086 s
accepted curved shock area    = 2.67211070324907966e20 m2
electron incident rate        = 3.25066729609523584e35 s-1
electron macroparticle weight = 6.58034390853486657e32
configuration fingerprint     = c5ffd8a154b2277f
```

The migrated decks are intended to resolve the same physical values, but
the configuration fingerprint is computed from the normalized option record,
which now carries contract level 6 (the one-dash mode used an internal
level 4, the 2026-10-10 generic-format decks schema 5) and the explicitly
written former defaults; the fingerprint
`c5ffd8a154b2277f` is therefore expected to change.  These numbers are not
evidence for the current files.  Native status of the migrated decks
(2026-10-10, schema-5 spellings, `sep3d` build with the single compiled
species `H_PLUS`; the native rerun with the schema-6 spellings is pending,
while `CFG3D18` confirms that all three run files parse):

* `amps.in` read as the no-argument default (`./amps --initialization-only`
  in a copy of this directory) completes initialization; products are in
  `test_output/input-migration/amps-in-default/`.
* `amps-injection-smoke.in` (rebuilt after the species-matching fix; logs in
  `test_output/srcsep3d-shock-injection/hplus-rank{1,4}/run.log`) gives the
  identical receipt on one and four ranks:

  ```text
  event fingerprint             = e82d70fa…29b3
  dt (global_time_step_s)       = 2.02430556841030818 s
  normalization radius / epoch  = 3.54807e10 m / 496.928571430180 s
  H_PLUS macroparticle weight   = 4.27644841713363414e32
  configuration fingerprint     = b882ccfb611c1b9a
  step 0: accepted_faces=1632 rate=2.12844622696998611e35 s-1 candidates=993  injected=993  disconnected=0
  step 1: accepted_faces=1632 rate=2.12838041717509552e35 s-1 candidates=1064 injected=1064 disconnected=0
  ```

  The compiled `H_PLUS` (AMPS mass 1 u) is identified as the proton
  population by charge number and mass number (`IdentifyPopulation` in
  `runtime/particle_normalization.cpp`, test `CFG3D14`); the CODATA-mass-only
  match that rejected it before is gone.  The smoke deck normalizes at the
  handoff radius, so its numbers are not comparable with the historical
  `ELECTRON` receipt below (which normalized at 10200 s).

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

Before the input migration the recorded results were respectively `PASS=151 FAIL=0 SKIP=0 ERROR=0`
and `PASS=57 FAIL=0 SKIP=0 ERROR=0`; the complete application aggregate is
`PASS=173 FAIL=0 SKIP=29 ERROR=0`. Its 29 SKIPs are declared external/native/
observational evidence gates, not failures of this source. Current native
one-/four-rank evidence is recorded in `CODEX_REDUCED_SHOCK_PLAN.md`; the
earlier initialization-only runs do not qualify native injection.  After the
migration, `CFG3D13` covers the unified input (includes, continuation,
derived-numerics keys, CLI spellings, diagnostics).  Recorded for the
migrated tree on 2026-10-10: `test-stage1` 153/153, `shock-front-test`
57/57, and `run_tests.py --all` `PASS=175 FAIL=0 SKIP=29 ERROR=0` (inside
`srcSEP3D/test/run_all_test.sh`).

Edit the `run` block of `amps.in` (or `amps-injection-smoke.in`) to change
the finite iteration horizon (`run.maximum_time_steps`), maximum represented
speed (`run.maximum_particle_speed_m_s`) or margin
(`run.time_step_margin_factor`).  Edit its `source: shock` block to change
the target mean model-particle count at the normalization radius
(`source.shock.samples_per_step`), kinetic-energy interval, phase-space power
law (`source.shock.spectrum_model`, `source.shock.phase_space_power_index`),
fatal Poisson guard (`source.shock.maximum_events_per_species_per_step`) or
required source-normalization radius (`source.shock.normalization_radius_m`).
Edit `sep3d.in` for the mesh, transport, turbulence, observer and output
settings shared by both runs.  Edit the `physics: shock_front` block in
`reduced-shock-surface.in` to change launch, handoff, ambient or magnetic
parameters.  In this derived mode `run.time_step_s`,
`source.shock.physical_particle_rate_per_s`,
`source.shock.injection_efficiency` and `species.macroparticle_weight` must
not be written; the parser rejects them. Any magnetic asset
change requires updating `assets.harmonics_sha256`; an incorrect checksum is a
fatal data-integrity error. The normalization radius must lie on the declared
monotonically outward trajectory and contain a positive accepted-shock area.

Startup reports `dt=f h_min/v_max` and, for each compiled species,
`W_s=Ndot_s dt/N_model`. The gross `Ndot_s` is the accepted upstream incident
flux and no hidden acceleration efficiency is applied. At each step the AMPS
user-defined injection callback recomputes accepted triangular-face rates and
draws exponential waiting times at `lambda=Ndot_live/W_s`. Therefore
`source.shock.samples_per_step` is the expected count only at the normalization
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
