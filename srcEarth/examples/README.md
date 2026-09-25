# Standalone Step 7 examples

`standalone_step7_combined.in` exercises the completed standalone product path with an
analytic dipole. From the AMPS repository root, run:

```bash
mpirun -np 4 ./amps -mode gridless \
  -i srcEarth/examples/standalone_step7_combined.in \
  -gridless-threads 16
```

The single invocation requests both product families from one immutable field snapshot.
It writes:

- `standalone_run_manifest.json` with model, field representation, snapshot identity,
  authoritative epoch, driver-validation status, and selected products;
- `cutoff_gridless_points.dat` plus
  `cutoff_gridless_dir_access_point_0000.dat` for cutoff and `A(E,Omega)`;
- `gridless_points_spectrum.dat` for boundary/local nominal and uncertainty spectra;
- `gridless_points_density.dat` for number density; and
- `gridless_points_flux.dat` for total, configured-channel, planar, and synthetic
  detector-response products.

The input enables the unresolved-trajectory failure gate. A nonzero exit must be
investigated; do not turn that gate off to make a production or validation run pass.
The modest energy/angular grids keep the example tractable and are not a convergence
claim. Use F/C validation profiles and resolution sweeps for science production.

To use a phenomenological field, replace the `#BACKGROUND_FIELD` block with one of the
released models and its complete inline parameters or a validated `DRIVER_FILE`:
T96, T01, T05/TS05, TA15N, TA15B, or TA16. File-backed tables must cover every requested
epoch, contain the exact model-specific columns, and use the documented native units.
The run fails before tracing if any condition is not met.

Mode3D uses the same input physics after changing `FIELD_EVAL_METHOD` in a copy of the
input to `GRID_3D`:

```bash
mpirun -np 4 ./amps -mode 3d \
  -i my_step7_combined_mode3d.in \
  -mode3d-threads 16
```

For meaningful gridless/Mode3D comparison, explicitly set a Mode3D mesh-resolution
profile and run the three-resolution I-F04 convergence study. Step 7 is standalone-only;
the SWMF-coupled product path is developed in Roadmap Steps 9–11.
