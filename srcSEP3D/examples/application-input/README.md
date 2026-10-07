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

The final 2026-10-07 clean build was exercised with one and four ranks. Both
layouts produced the same immutable receipt:

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
```

The recorded final results are respectively
`PASS=150 FAIL=0 SKIP=0 ERROR=0` and
`PASS=57 FAIL=0 SKIP=0 ERROR=0`.

Edit `sep3d.in` to change the model-particle count, maximum represented speed,
CFL margin or required source-normalization radius. Edit the reduced subsection
to change launch, handoff, ambient or magnetic parameters. Any magnetic asset
change requires updating `assets.harmonics_sha256`; an incorrect checksum is a
fatal data-integrity error. The normalization radius must lie on the declared
monotonically outward trajectory and contain a positive accepted-shock area.

Startup reports `dt=f h_min/v_max` and, for each compiled species,
`W_s=Ndot_s dt/N_model`. The gross `Ndot_s` is the accepted upstream incident
flux. It does not represent acceleration efficiency, and this example does not
activate reduced-front particle injection. The provider remains ambient-only
away from the surface: no output value should be interpreted as a sheath or
ejecta plasma state.
