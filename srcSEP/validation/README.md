# Step 15 scientific-validation campaign

Step 15 adds scientific verification after the transport architecture and
component contracts have stabilized. It deliberately reports five evidence
classes separately. Passing a manufactured solution is not called coupled
validation, and passing coupled code is not called observational validation.

## Included cases

| Case | Evidence class | Reference | Current status/prerequisite |
|---|---|---|---|
| `CV01` | Linked-application controlled numerical verification | Linked srcSEP/AMPS `--test CV01` output versus an independent closed-form ballistic characteristic with periodic/open boundary solutions | Implemented; requires linked executable |
| `VAL01` | Numerical verification | Independent analytical diffusion moments and adiabatic-cooling characteristic | Implemented |
| `VAL02` | Cross-mover verification | Matched `fte-dmumu`/`fte-mfp` mean-free-path closure | Implemented |
| `VAL03` | Cross-model verification | Independently coded conservative finite-volume pitch-angle solver | Implemented |
| `VAL04-SWCME` | Coupled integration | Actual `swcme::sep::Interface1D` states consumed by the srcSEP Parker kernel | Implemented |
| `VAL04-SWMF` | Coupled integration | Checksum-verified native SWMF replay | Requires external run evidence |
| Event campaigns | Observational validation | Held-out spacecraft products with uncertainties and forward operators | Requires external data evidence |

The source-only command is:

```sh
make test-scientific-validation
```

The new numbered validation portfolio uses the shared case registry described
in [cases/README.md](cases/README.md). CV01 is the first complete template:

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV01 \
  --output-dir /absolute/path/to/evidence/CV01
make test-cv01-unit SEP_EXECUTABLE=/absolute/path/to/amps
```

This verifies that the supplied executable advertises CV01, executes every
numerical stage through that linked application's native registry, and produces
the model input snapshot, independent reference, native JSON/JUnit, particle/
moment/profile CSV files, frozen metrics, negative-control evidence, PNG/EPS
figures, logs, and executable provenance. Later `CV02`–`EV02` cases are added
to the same registry and follow the same linked-application lifecycle.

It runs `VAL01`–`VAL04-SWCME` under AddressSanitizer and
UndefinedBehaviorSanitizer, assembles a temporary campaign, and proves that
release mode refuses the expected missing SWMF and spacecraft evidence. The
normal command leaves no generated result files in the source tree. To retain
its JSON, JUnit, and campaign reports outside the tree:

```sh
STEP15_REPORT_DIR=/absolute/path/to/evidence make test-scientific-validation
```

For a bounded equation-level regression that does not enter Step 15 campaign
assembly, use `make test-controlled-analytical`. Native execution is a separate
gate and requires the linked application:

```sh
make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps
```

For review figures, run the same controlled catalog through the Python
orchestrator. It preserves the component JSON and produces PNG plus EPS without
changing the scientific status or evidence class:

```sh
python3 test/run_tests.py --suite controlled-analytical \
  --output-dir /absolute/path/to/evidence/controlled-with-figures
```

Figures backed by a numerical/analytical CSV are pointwise solution overlays.
Where a case emits only moments, norms, convergence order, or conservation
residuals, the figure is explicitly a metric-versus-reference/acceptance plot.
The latter is useful for review but must not be cited as a spatial or temporal
solution curve. `analytical_plot_manifest.json` records the kind and source of
every figure; registry JSON remains the authoritative result.

WP11--WP20 also provide a bounded coefficient and stochastic-numerics gate:

```sh
make test-wp11-wp20-unit
```

That command checks bounded pitch diffusion, error controls, pure source-bound
coefficient kernels, species normalization, adaptive quadrature, and ballistic
compatibility. Its PASS is dependency-light numerical evidence only; it cannot
satisfy the native AMPS, real SWMF, or held-out observational gates below.

WP31--WP41 add a source-contract gate and a registry-generated compatibility
inventory:

```sh
make test-wp31-wp41-unit
make print-configuration-matrix
```

Evidence is labelled `analytical-core`, `source-integration`, `native-amps`,
`swmf-replay`, or `observational-validation`. A passing claim may not name a
level above what actually executed. The source gate therefore leaves native
adapter execution, scheduled multi-seed campaigns, MPI/OpenMP scaling, real
SWMF, and held-out observations explicitly BLOCKED when their inputs are absent.

## Completing external gates

Copy the appropriate file from `manifests/`, replace every placeholder, and
archive each input alongside the manifest or at a path relative to it. Every
input record must include its lowercase SHA-256, role, source URL or persistent
identifier, access time, and licensing/acknowledgment text. The runner opens
the local bytes and recomputes every checksum; a citation without archived
bytes is insufficient.

Each external class can be authenticated independently before assembling the
full release campaign:

```sh
make test-swmf-validation SWMF_MANIFEST=/evidence/swmf-replay.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS="/evidence/event-1.json /evidence/event-2.json"
```

Both targets require every manifest to declare `status=PASS` and verify every
referenced input checksum. The SWMF target accepts only
`COUPLED_INTEGRATION`/`SWMF_OUTPUT`. The observational target accepts only
`OBSERVATIONAL_VALIDATION`/`SPACECRAFT_OBSERVATION` and requires the supplied
held-out events, collectively, to cover all six required metric families.
These commands do not run or inherit status from controlled or native tests.

```sh
python3 validation/run_campaign.py \
  --numerical-json /evidence/step15-numerical-results.json \
  --swcme-json /evidence/step15-swcme-results.json \
  --swmf-manifest /evidence/swmf-replay.json \
  --observational-manifest /evidence/event-1.json \
  --observational-manifest /evidence/event-2.json \
  --output-dir /evidence/release --campaign-id campaign-identifier --release
```

`--release` succeeds only if all numerical cases pass, both the real SWCME and
real SWMF coupled cases pass, and the held-out observational manifests jointly
cover onset, anisotropy, spectra, fluence, decay, and multi-spacecraft
longitude. Any missing class is `INCOMPLETE`; an asserted failure is `FAILED`.

## Provenance and output contract

`campaign.schema.json` describes the assembled
`srcsep-validation-campaign-v1` record. `run_campaign.py` captures:

- the deterministic source-tree SHA-256 and every input/report SHA-256;
- complete internal case configuration, seeds, metrics, tolerances, and units;
- compiler command/version/flags and host platform;
- external run configuration, data provenance, output variables and units;
- independent statuses for numerical, cross-mover, cross-model, coupled, and
  observational evidence.

Both JSON and Markdown campaign reports are written transactionally through a
same-directory temporary file, flush, `fsync`, and atomic replacement. The JSON
record is authoritative; Markdown is a compact status view.

## Scientific-claim safeguards

An observational manifest must declare `data_class=SPACECRAFT_OBSERVATION`,
identify at least one mission/instrument data product, be marked held out,
provide an uncertainty method, and pass checksum verification. Labels that
identify manufactured data are rejected. This safeguard test is itself only a
software-governance test and never appears as observational evidence.

WP39 also requires `configuration.forward_operator` to name a versioned energy
response, angular response, cadence, species, dead time, saturation policy,
background subtraction, and uncertainty propagation. The C++ forward model
uses SI differential intensity, energy edges, `m² sr` geometric factor, and
seconds to predict detector counts. A synthetic response test validates that
operator, but only checksum-verified held-out spacecraft products can satisfy
the observational gate.

See each `cases/VAL*/README.md` for equations, configurations, metrics, and
limitations. See [STEP15_VALIDATION_REPORT.md](../STEP15_VALIDATION_REPORT.md)
for the evidence obtained in this delivered environment.
