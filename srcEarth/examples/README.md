# Standalone and coupled product examples

## Step 11 coupled flux and spectra

`swmf_step11_flux_spectrum.in.template` is the live coupled configuration for the
Step-11 transaction. It requests cutoff/access and density/spectrum products together,
uses the released static-magnetic characteristic, retains the 1% unresolved gate, and
declares example energy channels and a compact detector response. Replace every
`REPLACE_*` token; the file is intentionally not runnable with placeholders.

At each accepted PT epoch the coupled bridge selects matching trajectory/ephemeris
samples, collectively fingerprints the spectrum, species, energy grid, channels,
responses, and observation state, then calls the same `RunDensityAndFlux()` integrator
used by standalone Mode3D. A successful epoch writes close-verified density, local
spectrum, integral flux/channel/response, and termination files plus
`swmf_flux_spectrum_manifest<LIVE_SUFFIX>.json`.

Replay the exported field with `standalone_step9_swmf_replay.in.template`, copying the
complete live spectrum/time policy, ordered channels, response definitions, output
domain, and numerical controls. Compare all artifacts at once:

```bash
python3 srcEarth/test/USWMFCoupledProducts/compare_flux_products.py \
  --live-manifest coupled/swmf_flux_spectrum_manifest.swmf_t..._sid....json \
  --replay-dir replay \
  --json step11_live_replay_comparison.json
```

The default comparison is exact and rejects missing products, changed response or
observation schemas, different field/spectrum epochs, non-closing termination counts,
or unresolved fractions above the retained gate. The linked convergence, parallel,
instrument-closure, and O3 observation procedures are specified in
`../test/USWMFCoupledProducts/README.md`. The template remains quasi-static and
magnetic-only; it does not release an electric/time-dependent characteristic.

## Step 10 coupled cutoff/access

`swmf_step10_cutoff_access.in.template` is the live half of the Step-10
live/export/replay validation. It is consumed by an AMPS component launched inside
SWMF/PT, not by a standalone `./amps` command. Replace every `REPLACE_*` token with a
campaign value and archive the resolved input beside the field snapshot, status,
manifest, and cutoff/access files.

The callback cadence comes from `#TEMPORAL / FIELD_UPDATE_DT` in minutes, but the
scheduler compares the authoritative PT simulation clock in seconds. All ranks make
the same RUN/SKIP/DUPLICATE/STALE decision. A successful run writes a suffix containing
the exact nine-decimal simulation time and full content-derived snapshot ID. Do not
rename only part of an artifact set: the suffix is the transaction key shared by the
field export, cutoff/access output, diagnostic mesh file, status, and manifest.

Choose one explicit physical escape contract in `#DOMAIN_BOUNDARY`. `BOX` preserves
six-face escape. `SHUE` uses the Shue magnetopause plus the negative-X tail cap; the
configured Cartesian box still limits available mesh data, but a non-tail face inside
the magnetopause is a calculation failure rather than physical access. For `AUTO`
Shue coefficients, retain the exact PDYN and IMF Bz values used to resolve the surface.

After the live run, fill `standalone_step9_swmf_replay.in.template` with the exact
exported metadata and unchanged particle, output-domain, mover, search, and failure
controls. Compare the resulting table with:

```bash
python3 srcEarth/test/USWMFCoupledAccess/compare_cutoff_access.py \
  --live cutoff_3d_points<LIVE_SUFFIX>.dat \
  --replay cutoff_3d_points<REPLAY_SUFFIX>.dat \
  --json step10_live_replay_comparison.json
```

The default comparison is exact and first requires matching epoch, snapshot ID,
content fingerprint, mesh revision, and boundary policy. Nonzero tolerances must be
explicitly justified by the registered I-F04/I-F05 study; they are not a way to bypass
a provenance mismatch. Repeat the complete live/replay case at 1x1, 2x8, and default
8x16 and after restart at the same epoch. This release is quasi-static and
magnetic-only; it does not validate a time-dependent electric characteristic.

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
