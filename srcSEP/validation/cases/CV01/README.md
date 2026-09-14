# CV01 — ballistic streaming on a uniform field line

## Why this case is required

CV01 isolates the irreducible characteristic of field-line SEP transport. A
failure cannot be attributed to turbulence, diffusion, focusing, shock
injection, cooling, or background interpolation because those operators are
explicitly disabled. The test therefore detects a reversed pitch-angle or
field-line convention, an incorrect relativistic speed, an erroneous streaming
normalization, timestep-dependent packet broadening, biased boundary timing,
or lost statistical weight.

CV01 is a linked-application numerical-verification gate. Every numerical
realization is executed by the requested srcSEP/AMPS binary through its native
`--test CV01` registry entry, which calls the compiled-in
`AdvanceFocusedTransportDmumu()` production core. CV01 deliberately exits
before AMPS mesh initialization because its uniform synthetic line must remain
independent of production background setup; it validates executable linkage,
production CLI/registry dispatch, the compiled transport kernel, boundary
accounting, output, and the external comparison workflow. It does not claim
SWMF coupling or observational validation.

## Governing solution

For a straight line, constant plasma advection `U`, uniform `|B|`, zero
`D_mumu`, and no velocity gradients, each particle satisfies

```text
s(t) = s0 + (U + mu*v) t
mu(t) = mu0
p(t) = p0
```

The independent Python reference converts the configured proton kinetic energy
to speed and momentum relativistically, evaluates this characteristic, and
applies either a periodic coordinate map or the exact first crossing of the
open interval. It never imports or calls srcSEP numerical code.

## Reviewed default input

`input.json` is the complete model input. It specifies SI constants, a 1 AU
uniform line, 5 nT magnetic field, 400 km/s parallel solar wind, 5 cm^-3 number
density, 10 MeV protons, pitch-angle cosines `-1`, `-0.5`, `0.5`, and `1`, and a
weighted Gaussian packet centered at 0.5 AU. Periodic and open boundaries are
both run for 1800 s with 60, 30, and 15 s timesteps and a 60 s sampling cadence.

Every non-streaming operator is recorded as disabled. The input validator
rejects a nonzero magnetic gradient, flow gradient/divergence, `D_mumu`, source,
loss, perpendicular transport, or particle-wave coupling instead of silently
changing the meaning of CV01.

## Execution

Run CV01 through the common user-facing runner:

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV01 \
  --output-dir test_output/CV01
```

The path may instead be supplied once through `SEP_EXECUTABLE`:

```sh
SEP_EXECUTABLE=/absolute/path/to/amps \
  python3 test/run_tests.py --validation-case CV01 \
  --output-dir test_output/CV01
```

To validate an AddressSanitizer/UndefinedBehaviorSanitizer build, build the
linked application with those flags and pass that binary through `--amps`.
The runner never recompiles or substitutes a model driver.

Run a reviewed input variant without changing the catalog default:

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV01 \
  --case-input /absolute/path/to/CV01-input.json \
  --output-dir test_output/CV01-variant
```

The bounded target always checks that the compiled-in CV01 source builds with
strict warnings. It runs the complete application validation only when
`SEP_EXECUTABLE` names the linked binary; otherwise it reports the native gate
as `SKIP`:

```sh
make test-cv01-unit SEP_EXECUTABLE=/absolute/path/to/amps
```

`validation/run_case.py --amps /path/to/amps --case CV01` is the lower-level registry runner. It is
useful for automation that wants raw reports without the common overlay plot;
normal users should prefer `test/run_tests.py`.

## Model, reference, and negative-control separation

- `cv01_model.cpp` is compiled into the linked application. Its native registry
  callback reads the immutable packet fixture and advances every state with the
  production focused-transport core. A provider satisfying the full
  coefficient/provenance contract returns exactly zero diffusion. Its utility
  headers are referenced relative to the source file so the normal AMPS build
  does not require an extra `srcSEP/util` include-directory modification.
- `reference_solution.py` evaluates the closed-form characteristic in a
  separate Python process. It can also be invoked directly with `--help`.
- `case.py` generates the shared initial packet, runs six linked-application and
  reference pairs (two boundaries times three timesteps), verifies every native
  JSON report, reloads saved CSVs, computes moments and acceptance metrics, and
  renders the prescribed four-panel plot.
- The negative control deliberately compares the `mu=1` result with a
  sign-reversed reference. CV01 fails if this wrong solution is not separated
  from the model by at least the configured fraction of the line length.

## Acceptance metrics

The default gate requires maximum particle and centroid error below `1e-12` of
the line length; variance error below `1e-12` of line length squared; relative
weight closure below `1e-13`; crossing-time error below `1e-12` of final time;
and pitch-angle and relative-momentum errors below `1e-13`. All four pitch
angles, both boundaries, all saved times, and all three timesteps contribute to
these maxima.

Constant-coefficient ballistic streaming is integrated exactly by the
production midpoint characteristic, so a finite convergence order is not
mathematically measurable once roundoff dominates. The three timesteps instead
form an exactness/refinement-invariance gate and expose any end-of-step boundary
time bias. Nontrivial temporal order is measured by later manufactured cases.

## Output contract

The case output directory contains:

- `resolved_input.json` and `initial_particles.csv`;
- `native/<boundary>_dt_<dt>/CV01_model.csv`, native JSON/JUnit, and generated
  native argument manifests for every linked realization;
- `reference_*.csv` for every independent boundary/timestep solution;
- `CV01_solution.csv`, the standard centroid overlay series;
- `CV01_packet_moments.csv` and `CV01_final_profiles.csv`;
- `CV01_four_panel.png` and `.eps`, showing representative packet profiles,
  centroid error, variance error, and retained/escaped statistical weight;
- `CV01_run.log` and `provenance.json` with exact linked commands and SHA-256
  checksums of the executable, source, configuration, and outputs.

The run fails before generating a reference if `--amps` is missing, is not
executable, rejects `--list-tests`, or does not advertise `CV01`. It also fails
if any native process exits nonzero, its registry result is not `PASS`, or its
model CSV is absent. This prevents a Python-only result from being confused
with validation of the linked application.

A generic native `--all-tests`, `--test-group controlled-analytical`, or direct
`--test CV01` command without the case paths reports CV01 as `SKIP`; it does not
invent defaults or make unrelated registry runs fail. Use the Python validation
command above for the complete case.

The parent directory also contains `srcsep-tests.json`,
`validation-run-manifest.json`, `run_manifest.json`, and the common
`plots/CV01_comparison.png/.eps` overlay with a residual panel. JSON metrics are
authoritative; figures are review aids generated from saved CSV artifacts.

## Failure interpretation

- A signed trajectory error that reverses with `mu` indicates field-line or
  pitch-angle orientation.
- An energy-dependent displacement error implicates relativistic conversion.
- Packet variance growth indicates unintended diffusion or inconsistent
  coordinate/bin handling.
- Timestep-dependent crossing error indicates end-of-step rather than exact
  boundary accounting.
- Active plus escaped weight differing from unity indicates deletion/sampling
  bookkeeping loss.
- A failed negative control means the test lacks power and must not be accepted
  even when the nominal model/reference comparison is small.
