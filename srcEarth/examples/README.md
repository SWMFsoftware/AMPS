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
profile and run the three-resolution I-F04 convergence study.

## Step 9 SWMF snapshot replay

`standalone_step9_swmf_replay.in.template` is the offline half of the Step-9
live/replay comparison.  First run the coupled case with
`SWMF_SNAPSHOT_EXPORT T`, `SWMF_SNAPSHOT_EXPORT_PREFIX <stem>`, and the default
`SWMF_DERIVED_ELECTRIC_FIELD OFF`.  Preserve the emitted CSV, product status JSON,
input, executable revision, and all live product files together.

Copy the template and replace every `REPLACE_*` token with the corresponding live
value.  In particular, reproduce the exact absolute epoch, Cartesian domain, AMR
resolution controls, product domain, particle/spectrum definition, energy and angular
grids, mover budgets, and unchanged failure thresholds.  `SWMF_SNAPSHOT_FILE` must
name the exported CSV.  Do not use gridless mode: replay is intentionally an exact
Mode3D mesh-state comparison rather than refitting the SWMF state to an analytic field.

Run the replay from the repository root, for example:

```bash
mpirun -np 4 ./amps -mode 3d \
  -i srcEarth/examples/my_step9_swmf_replay.in \
  -mode3d-threads 16
```

Startup fails before tracing if the constructed mesh differs in block dimensions,
leaf count, cell keys, centres, or mesh revision, or if schema, GSM/SI units, domain,
epoch, electric mode, content fingerprint, or snapshot ID is inconsistent.  Acceptance
requires pointwise comparison of cutoff, directional access, differential spectra, and
integral flux against the live products at the existing plan thresholds, then a repeat
after restart and across 1x1, 2x8, and default 8x16 rank/thread layouts.  Missing output
or a non-PASS status is a failure; it is never an excluded observation.

## Step 8 event-campaign example

`standalone_step8_campaign/` packages the Step-7 analytic dipole case as a frozen,
restartable three-epoch workflow. Its manifest includes every required resource class,
provenance statement, SHA-256 digest, exact units, expected Tecplot variable order and
row count, and the unchanged unresolved-support gate. The synthetic identical
east/west responses and unity observation are an analytic orchestration reference;
they are not a production O1/O2 validation result.

Check the complete package without an AMPS executable or network access:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --output-dir test_output/step8_preflight \
  --validate-only
```

Render all epoch inputs without executing them:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --amps ./amps \
  --output-dir test_output/step8_dry_run \
  --dry-run
```

For a real run, replace `--dry-run` with `--restart`. A previous run is skipped only
when its fingerprint and every recorded artifact hash still match. Production event
manifests must replace the synthetic resource files with independently archived O1,
O2, or preregistered holdout resources while retaining the same fail-closed schema and
validation gates.
