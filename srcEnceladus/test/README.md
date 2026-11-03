# Enceladus forward-model tests

Run the complete application-independent suite from `srcEnceladus`:

```bash
python3 test/run_tests.py
```

Run one case, with its own JSON/SVG receipt, by its stable ID:

```bash
python3 test/run_case.py U06_source_law
python3 test/run_case.py U08_validation
```

`python3 test/run_case.py --help` lists all exporter, compile, parser,
chemistry, source-law, numerics, and negative-validation cases.

Select the application through the repository workflow (never by editing
`.amps.conf`), confirm no Mars2 AMPS build/run is active, remove only the
generated `build` directory, regenerate, and compile before validating the
pre-MPI startup gates:

```bash
cd /home/vtenishe/Mars2/AMPS
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_native_startup.py \
  --executable /home/vtenishe/Mars2/AMPS/amps
```

This checks CLI help, validation of the complete checked-in E3 deck, explicit rejection of a
deliberately corrupted inward copy of the projected stripe receipt, pre-MPI
rejection of an inaccessible configured plasma-cache directory using the
production projection, and final generated pointer-mover and loss-callback
macros. These tests do not initialize MPI, allocate a mesh, or claim an
evolving physics run.

The global runner executes the standalone source exporter tests, exports a
complete 8-jet/7-stripe fixture, builds the dependency-free C++ core, and runs
focused source-parser, runtime-parser, Burger-chemistry, Tenishev flux-law,
canonical ID/geometry, radial stripe projection, parent-loss-only planning,
Courant/weight, conditional flyby-mesh-refinement, strict inactive-control, and dimensionless
velocity-shift population-control tests. The source-law case compares the production K3
moment and inverse CDF with direct composite-Simpson quadrature; it does not
reuse the production recurrence. The chemistry case verifies that Burger's
resonant channel contains both H2O and H2O+. It then compares named outputs with
`reference/forward_model_reference.json` using independent equations in the
Python runner. Production code never writes the expected values.

Results are `PASS`, `FAIL`, or an explicit compilation/runner failure in
`test_output/srcEnceladus/test_receipt.json`. The adjacent
`reference_comparison.svg` plots relative errors and is generated with the
Python standard library, so visualization does not depend on matplotlib.

The fixture records a 252.100 km minimum and maximum radius for every exported
stripe node and the same explicit reference radius. The exporter unit test
independently recomputes every vertex norm; the native negative case proves
that a one-percent inward corruption is still rejected.

The exporter tests also validate both E3 and E5 H2O+CO2 startup manifests:
each species has all 15 equal strengths, the sums reproduce the selected H2O
and apparent-CO2 totals exactly, and the CO2 provenance retains the factor-2--5
fragmentation bracket and its derived physical-rate interval.

Run the bounded one-rank E3 integration and conservation gate in a new output
directory (the runner refuses to overwrite prior evidence):

```bash
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --output test_output/srcEnceladus/e3_first_run_qualified
```

The runner requires normal process completion, the exact E3/no-plasma receipt,
a positive global timestep, 800 completed steps and injection calls, exactly
80,000 injected particles of each primary, zero CO2 photolytic losses, no
plasma-output file, native finalization, and exact per-species conservation:

```text
injected = photolytic loss + surface loss + outer escape + resident
```

`--launcher direct` is available only for one-rank MPI singleton diagnostics.
Some restricted containers deny the local sockets needed by both launcher
forms; that is an MPI-environment failure before model startup, not permission
to skip the native run. Use an authorized host MPI environment and retain the
failed environmental receipt separately.

The 2026-11-01 qualified E3 receipt passes all nine gates, but this does not
qualify four-rank ownership, plasma/coupled chemistry, numerical convergence,
or observational accuracy. Those integration gates remain in
`development_and_validation_plan.md`.

For the refined 10,000-step run, the runner additionally requires
`flyby_mesh_refinement=on` in the native receipt and verifies
`cell_max_edge_m <= 20000` for every in-domain row in every Tecplot zone. This
is a direct check of the allocated mesh at the observation locations; it does
not infer resolution from the requested configuration.

## Repeatable INMS and UVIS comparison run

The checked-in runtime deck requests the currently usable INMS flybys and UVIS
occultations through repeatable named `observation` blocks. E7 supplies a
calibrated apparent-CO2 density comparison; E13 remains model-only; E21
retains H2O/CO2 detector counts as display-only native context because no
validated scalar density response exists. Additional UVIS events without
qualified reductions/rays are deliberately not treated as comparison data.
To add a new data set, add another block; to omit one from a run, comment its complete
`#begin observation` through `#end observation` range. The parser contract,
units, required columns, ray geometry, and fail-closed validation rules are
documented in `../RUNTIME_INPUT.md`.

Regenerate the repository-released E7/E21 tables and their checksummed
provenance manifest from the validation archive with:

```bash
python3 standalone/export_forward_observations.py
```

Run the current E3 source distribution for ten 1,000-step sampling windows in
a fresh output directory:

```bash
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --input /home/vtenishe/Mars2/AMPS/input/enceladus.in \
  --output test_output/srcEnceladus/e3_10000_extended_observations \
  --launcher direct --timeout-s 10800 \
  --expected-iterations 10000 --expected-sample-length 1000 \
  --expect-volume-fields off
```

Each configured request produces a ten-zone Tecplot file named
`PT/plots/observation.<request>.dat`. Generate the final-window comparison
PNGs and the machine-readable summary with:

```bash
python3 test/plot_observation_comparisons.py \
  test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/observation.*.dat \
  --output-dir test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/comparisons \
  --summary test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/comparisons/observation_comparison_summary.json \
  --expected-zones 10
```

The plotter overlays observations only when the input block declares a
response-qualified `comparison_species`. A request with
`comparison_species=none` is either model-only or, when the three
`context_*` keys are present, adds a third native-observable panel. Context is
never scored against model density and is labelled explicitly in the plot and
JSON receipt.
UVIS panels compare the simulated line-of-sight H2O column with the optical-
depth shape using one fitted nonnegative scale plus a signed baseline, and
separately report the AMPS peak/published-column ratio. A plotter `PASS` means
that all products have the requested zone count, stable row sequence, valid
schema, and a generated PNG; it is not a claim of observational agreement or
physical convergence.
