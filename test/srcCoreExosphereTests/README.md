# Native spherical exosphere-source qualification

This application exercises the actual legacy
`Exosphere::SourceProcesses::InjectionBoundaryModel` on a one-rank 3-D AMR
mesh.  It is intentionally separate from dependency-light source-law tests:
the latter validate exact inverse CDFs and normalization, while this native
receipt validates particle allocation and AMR/surface integration.

The generated profile compiles the five spherical modes that can be tested
without a target-specific source adapter: impact vaporization, photon
stimulated desorption (PSD), thermal desorption, uniform-rate solar-wind
sputtering, and vertical injection.  The source table is checked at compile
time and every mode must emit at least one particle in the fixed random
stream.  External-domain and background-plasma injection are domain-boundary
modes, not spherical modes.  Uniform-local-temperature injection is omitted
because its generated macro requires an application temperature-function
receipt; user-defined sources are application adapters.  Those exclusions are
reported as explicit scope limits, not as passes.

The case checks:

- native particle-buffer allocation and exactly one compiled source ID per
  particle;
- a finite one-rank AMR owner and membership in the committed owner-cell list;
- the configured individual statistical-weight correction;
- closure among `CalculatedSourceRate`, the per-element sampled source-rate
  ledger, and `PIC::BC::ParticleProductionRate`;
- desorption only for PSD, thermal, and sputtering, whose generated source
  table marks them inventory-depleting;
- deterministic accommodation accounting (`return=W`, `sticking=W`,
  `adsorption=0.75W`); and
- the time-dependent surface update
  `N_after=N_before+adsorption-desorption`, followed by clearing the physical
  interval flux buffers.

Source occurrence counts are stochastic observability checks, not estimates
of the physical rates.  Exact conservation uses the represented weight of
the particles that were actually emitted.

## Exclusive build and run

The suite runner owns the generated AMPS configuration and therefore must run
exclusively.  From the repository root:

```sh
python3 src/test/run_tests.py --case CORE-EXOSPHERE-SPHERE-INTEGRATION-001
```

For a manual qualification, first verify that no build/test process is using
this checkout, then follow the repository preflight literally:

```sh
rm -rf -- build
./Config.pl -application=test/core-exosphere-sphere
env MAKEFLAGS="-j16" make amps
mpiexec -n 1 ./amps
```

The application is qualified only for one MPI rank.  It does not validate
multi-rank source statistics, target-specific thermal maps, plasma/yield-mode
sputtering, external-boundary injection, or non-spherical/resolved surfaces.
