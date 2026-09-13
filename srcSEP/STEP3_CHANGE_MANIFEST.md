# Step 3 change manifest

Step 3 corrects shared flux-tube geometry, unit handling, and source
normalization. It builds on the Step 1 CLI/test registry and Step 2 immutable
background snapshot and simulation clock.

## New files

- `FLUX_TUBE_GEOMETRY.md`
- `flux_tube_geometry.cpp`
- `util/sep_physical_units.h`
- `util/sep_flux_tube_geometry_core.h`
- `util/sep_flux_tube_geometry_core.cpp`
- `test/run_step3_tests.sh`
- `test/step3/test_flux_tube_geometry.cpp`
- `STEP3_CHANGE_MANIFEST.md`

## Principal modified files

- `sep.h`, `shock_injection.cpp`, `field_line.cpp`, and `main_lib.cpp`;
- `shock_analytical_model.cpp` and `solar_wind.cpp`;
- wave-energy, turbulence, growth, mover, sampling, and main driver sites that
  consume segment volume or boundary area;
- `makefile`, `README.md`, and `test/README.md`.

## Quick verification

```sh
make test-geometry-source-unit
```

No generated object, archive, sanitizer, coverage, or test-output file belongs
in version control.
